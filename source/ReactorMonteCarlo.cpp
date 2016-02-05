/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactorMonteCarlo.cpp
 * Author: chris
 * 
 * Created on December 17, 2015, 7:09 PM
 */

#include "ReactorMonteCarlo.h"
#include "MicroGeometry.h"
#include <sstream>
#include "InputDataFunctions.h"
#include <chrono>
#include <thread>
#include "InfiniteCompositeReactor.h"


ReactorMonteCarlo::ReactorMonteCarlo() {}

ReactorMonteCarlo::ReactorMonteCarlo(InfiniteCompositeReactor* reactor,const Real &starting_k_eff, const  std::string &run_directory)
{
    Real k_eff;
    Real prompt_removal_lifetime;
    Real k_eff_sigma;
    Real prompt_removal_lifetime_sigma;
    
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _reactor = reactor;
    
    getRawCriticalityParameters(k_eff,k_eff_sigma,prompt_removal_lifetime,prompt_removal_lifetime_sigma);
    _virtual_k_eff_multiplier = starting_k_eff / k_eff;  
    _current_k_eff = starting_k_eff;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = k_eff_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
    
    //Parameters that will come from the Monte Carlo Simulation
    std::pair<FissionableIsotope,Real> U235 = { FissionableIsotope::U235, 1 };
    _fission_tally_listing = { U235 };
}

/**
 * Adjusts the k_effective to match the desired starting reactivity
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::updateAdjustedCriticalityParameters()
{
    Real raw_k_effective, raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma;
    
    getRawCriticalityParameters(raw_k_effective,raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma);    
    
    _current_k_eff = raw_k_effective * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = raw_k_effective_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
}

/**
 * Get the raw k_effective from MCNP
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::getRawCriticalityParameters( Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma)
{

    //clean up the previous MCNP data files
    std::string clean_mcnp_command = "rm " + this->_run_directory + "mcnp-composite-fue[l-z].out " + this->_run_directory + "runtp[e-z] " + this->_run_directory + "srct[p-z]";
    exec(clean_mcnp_command);
    
    //create a symbolic link to the Doppler broadened cross sections
    std::string symbolic_link_command = "cd " + this->_run_directory + "; ln -s ../../../doppler-broadened-cs/otf*txt .";
    exec(symbolic_link_command);
    
    //create the MCNP input file
    std::string file_root = "mcnp-composite-fuel";
    std::string input_file_name = file_root + ".inp";
    std::string output_file_name =  file_root + ".out";
    this->createMCNPOutputFile(input_file_name);
   
    std::string command_line_log_file = "mcnp_run_log.txt";
    
    //Run the file
    #ifdef LAPTOP
    
    //We are just running MPI here
    std::string mcnp_path = "/media/chris/DoubleSpace/MCNP/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    std::string command = "cd " + this->_run_directory + "; mpirun -np  7 " + mcnp_path + " i=" + input_file_name + " o=" + output_file_name + " | tee " + command_line_log_file;
    exec(command);
    
    #elif PRACTICE_CLUSTER
    
    std::string mcnp_path = "/share/apps/mcnp/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    
    //Delete the old log output by overwriting it with a blank string (we are using the log to check to see if the qsub is done)
    std::string remove_command = "cd " + this->_run_directory + "; echo \"\" > " + command_line_log_file;
    exec(remove_command);

    //Run Submission Script with the created MCNP file
    std::string qsub_command = "cd " + this->_run_directory + ";qsub -pe orte 32 ../../../composite-fuel-submission-script.sh " + input_file_name + " " + output_file_name + " " + command_line_log_file;
    exec(qsub_command);
    //Constantly read the output file until it says mcrun done 
    std::string search_lock = "cd " + this->_run_directory + ";cat " + command_line_log_file + " | grep \"mcrun  is done\"";
    std::string is_done;
    
    do
    {   
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        is_done = exec(search_lock);
    }while(is_done == "");
    
    #endif
    

    //Read the output file
    this->readOutputFile(output_file_name, k_eff, k_eff_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma);
        
}

void ReactorMonteCarlo::readOutputFile(const std::string &file_name, Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma)
{
    
    std::string base_command = "cd " + this->_run_directory + ";cat " + file_name + " | perl -ne ";   
    
    std::string k_eff_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $1'";
    std::string k_eff_command = base_command + k_eff_regex;
    std::string k_eff_str = exec(k_eff_command);
    k_eff = std::stod(k_eff_str);
        
    std::string k_eff_sigma_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $2'";
    std::string k_eff_sigma_command = base_command + k_eff_sigma_regex;
    std::string k_eff_sigma_str = exec(k_eff_sigma_command);
    k_eff_sigma = std::stod(k_eff_sigma_str);
    
    std::string prompt_lifetime_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $1'";
    std::string prompt_lifetime_command = base_command + prompt_lifetime_regex;
    std::string prompt_removal_lifetime_str = exec(prompt_lifetime_command);    
    prompt_removal_lifetime = std::stod(prompt_removal_lifetime_str);
        
    std::string prompt_lifetime_sigma_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $2'";
    std::string prompt_lifetime_sigma_command = base_command + prompt_lifetime_sigma_regex;
    std::string prompt_removal_lifetime_sigma_str = exec(prompt_lifetime_sigma_command);
    prompt_removal_lifetime_sigma = std::stod(prompt_removal_lifetime_sigma_str);
    
    
    /*std::ifstream mcnp_output_file;
    mcnp_output_file.open(file_name);
    
    std::string output_file_text;
    
    if (mcnp_output_file.is_open()) 
    {
        output_file_text.assign( std::istreambuf_iterator<char>(mcnp_output_file) , std::istreambuf_iterator<char>()  );
    }
    mcnp_output_file.close();
    
    
   
    //"estimated combined collision/absorption/track-length keff = 1.79363 with an estimated standard deviation of 0.00064"
    //the final combined (col/abs/tl) prompt removal lifetime = 2.6650E-04 seconds with an estimated standard deviation of 4.6932E-07 */
}

std::string ReactorMonteCarlo::getMaterialCards()
{
    std::stringstream material_cards, otfdb_card;
        
    Real enrichment_fraction = _reactor->_input_file_reader->getInputFileParameter("Uranium Enrichment Fraction", 0.2);

    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    size_t number_zones = geometry_data.size();   
       
    for( size_t index = 0; index < number_zones  ; index++ )
    {
        size_t current_zone = index + 1;
    
        Materials material = geometry_data[index].first; 
        std::string material_card_entry;
        std::string doppler_card_entry;
        _reactor->_micro_sphere_geometry->_material_library.getMcnpMaterialCard(material,current_zone,material_card_entry, doppler_card_entry, enrichment_fraction);
        
        material_cards << material_card_entry;
        otfdb_card << doppler_card_entry;        
    }       
    
     #ifdef LAPTOP

    std::string U238_cs = "92238.66c";
    std::string U235_cs = "92235.66c";

    #elif  PRACTICE_CLUSTER 

    std::string U238_cs = "92238.80c";
    std::string U235_cs = "92235.80c";    

    #endif
    
    otfdb_card << " OTFDB " << U238_cs << std::endl;
    otfdb_card << "       " << U235_cs << std::endl;
    //otfdb_card << "       8016.60c" << std::endl;
    //otfdb_card << "       6000.60c" << std::endl;
    
    
    return material_cards.str() + otfdb_card.str();
}

std::string ReactorMonteCarlo::getCellCards()
{
    std::stringstream cell_cards;
    
    //For now we will assume that the density stays constant regardless of temperature change to perserve conservation of mass
    Real density_derived_temperature = 400;
    const Real MeVperK = 8.617e-11;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    size_t number_zones = geometry_data.size();   
       
    for( size_t index = 0; index < number_zones  ; index++ )
    {
        size_t current_zone = index + 1;
    
        
        Real temperature = _reactor->_thermal_solver->getAverageTemperature(index);
        Materials material = geometry_data[index].first;
        //First is the density, second is the density derivative which we don't need. Divide by 1000 to convert from kg/m^3 to g/cm^3
        Real density = _reactor->_micro_sphere_geometry->_material_library.getDensityPair(material,density_derived_temperature,0).first/1000;
        
        
        
        //we want to put our reflecting boundary condition here
        if( current_zone == 1 )
        {
            cell_cards << " " << current_zone << " " << current_zone << " -" << density << " -"  << current_zone << " imp:n=1 TMP=" << temperature*MeVperK  << " $fuel T = " << temperature << " K" << std::endl;
        }
        else
        {
            cell_cards << " " << current_zone << " " << current_zone << " -" << density << " " << ( current_zone -1 )  << " " << " -"  << current_zone << " imp:n=1 TMP=" << temperature*MeVperK <<" $matrix T = " << temperature << " K" << std::endl;
        }
    }
    
    return cell_cards.str();
}

std::string ReactorMonteCarlo::getSurfaceCards()
{
    std::stringstream surface_cards;
    
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    size_t number_zones = geometry_data.size();   
       
    for( size_t index = 0; index < number_zones  ; index++ )
    {
        size_t current_zone = index + 1;
    
        //we want to put our reflecting boundary condition here
        if( current_zone == number_zones)
        {
            surface_cards << "*";
        }
        else
        {
            surface_cards << " ";
        }
        //cm sphere radius
        Real shell_radius = geometry_data[index].second * 100; 
        
        surface_cards << current_zone << " SPH 0 0 0 " << shell_radius << std::endl;
        
    }
    
    return surface_cards.str();
}


void ReactorMonteCarlo::createMCNPOutputFile(const std::string &file_name)
{
    
    std::stringstream mcnp_file;
   
     
    std::string cell_cards = this->getCellCards();
    std::string surface_cards = this->getSurfaceCards();
    std::string material_cards = this->getMaterialCards();
    int number_of_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles", 33);
    
    mcnp_file << "Composite Fuel Kernel Scale" << std::endl;
    mcnp_file << "c Simulating a small UO2 fuel kernel inside a graphite matrix" << std::endl;
    mcnp_file << "c TMP [MeV] = 8.617e-11 [MeV/K] * T [K]   " << std::endl;
    mcnp_file << "c ----------------------CELL CARDS----------------------------" << std::endl;
    mcnp_file << cell_cards;
    mcnp_file << "c end cell cards" << std::endl;
    mcnp_file << std::endl;
    mcnp_file << "c ----------------------SURFACE CARDS---------------------" << std::endl;
    mcnp_file << surface_cards;
    mcnp_file << "c end surface cards" << std::endl;
    mcnp_file << std::endl;
    mcnp_file << "c ------------------MATERIAL AND DATA CARDS--------------------" << std::endl;
    mcnp_file << material_cards;
    mcnp_file << " KCODE 10000 1.5 3 " << number_of_cycles << "  $need at least 30 active cycles to print results" << std::endl;
    mcnp_file << " KSRC 0 0 0" << std::endl;
    mcnp_file << " print" << std::endl;
    mcnp_file << "c end data" << std::endl;
    
    std::ofstream mcnp_input_file;
    mcnp_input_file.open( this->_run_directory + file_name,std::ios::out);
    mcnp_input_file << mcnp_file.str();
    mcnp_input_file.close();
    
    std::cout<<mcnp_file.str();
}
