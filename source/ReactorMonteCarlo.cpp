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
#include <math.h>
#include "InputDataFunctions.h"
#include <chrono>
#include <thread>
#include "InfiniteCompositeReactor.h"


ReactorMonteCarlo::ReactorMonteCarlo() {}

ReactorMonteCarlo::ReactorMonteCarlo(InfiniteCompositeReactor* reactor,const Real &starting_k_eff, const  std::string &run_directory)
{
    Real k_eff, prompt_removal_lifetime, k_eff_sigma, prompt_removal_lifetime_sigma;
    Real nd_k_eff, nd_prompt_removal_lifetime, nd_k_eff_sigma, nd_prompt_removal_lifetime_sigma;
    
    
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _reactor = reactor;
    
    //Cells per zone splits each zone up into multiple cells
    _cells_per_zone = this->_reactor->_input_file_reader->getInputFileParameter("Cells Per Zone", 1 );  
    
    

    getRawCriticalityParameters( "keffective-calc",            k_eff,    k_eff_sigma,    prompt_removal_lifetime,    prompt_removal_lifetime_sigma);
    getRawCriticalityParameters( "keffective-no-delayed-calc", nd_k_eff, nd_k_eff_sigma, nd_prompt_removal_lifetime, nd_prompt_removal_lifetime_sigma, false);    
    
    
    _virtual_k_eff_multiplier = starting_k_eff / k_eff;  
    _current_k_eff = starting_k_eff;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = k_eff_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
    _current_beta_eff = ( k_eff - nd_k_eff ) / k_eff;
    _current_beta_eff_sigma = this->getBetaEffSigma(k_eff, k_eff_sigma, nd_k_eff, nd_k_eff_sigma);
    

     
    
}

/**
 * Adjusts the k_effective to match the desired starting reactivity
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::updateAdjustedCriticalityParameters(const bool &update_beta)
{
    Real raw_k_effective, raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma;
    
    getRawCriticalityParameters("keffective-calc", raw_k_effective,raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma);    
    
    _current_k_eff = raw_k_effective * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = raw_k_effective_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
    
    if( update_beta )
    {
        Real raw_nd_k_effective, raw_nd_k_effective_sigma, nd_prompt_removal_lifetime, nd_prompt_removal_lifetime_sigma;

        getRawCriticalityParameters("keffective-no-delayed-calc", raw_nd_k_effective,raw_nd_k_effective_sigma, nd_prompt_removal_lifetime, nd_prompt_removal_lifetime_sigma);    
        
        _current_beta_eff = ( raw_k_effective - raw_nd_k_effective ) / raw_k_effective;
        _current_beta_eff_sigma = getBetaEffSigma(raw_k_effective, raw_k_effective_sigma, raw_nd_k_effective, raw_nd_k_effective_sigma );
    }
    
    
}

Real ReactorMonteCarlo::getBetaEffSigma(const Real &k_eff,const Real &k_eff_sigma,const Real &nd_k_eff,const Real &nd_k_eff_sigma)
{
    Real current_beta_eff = ( k_eff - nd_k_eff ) / k_eff;
    
    //sqrt of the sum of the squares of sigmas
    Real difference_uncertainty = sqrt( k_eff_sigma * k_eff_sigma + nd_k_eff_sigma * nd_k_eff_sigma );
    
    //Find the uncertainty fractions
    Real numerator_uncertainty_fraction = difference_uncertainty / ( k_eff - nd_k_eff );
    Real denominator_uncertainty_fraction = k_eff_sigma / k_eff;
    
    //uncertainty is the sum of the squares of the fractions
    return current_beta_eff * sqrt( numerator_uncertainty_fraction * numerator_uncertainty_fraction  + denominator_uncertainty_fraction * denominator_uncertainty_fraction );
}


/**
 * Get the raw k_effective from MCNP
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::getRawCriticalityParameters(const std::string &file_root, Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma,const bool &delayed_neutrons)
{
    //create the MCNP input file
    std::string run_title = this->_reactor->_run_name + " Time = "  + std::to_string(this->_reactor->_thermal_solver->_current_time);
    std::string input_file_name = file_root + ".inp";
    std::string output_file_name =  file_root + ".out";
    std::string runtpe_name =  file_root + ".runtpe";
    std::string srctpe_name =  file_root + ".srctp";
    
    //clean up the previous MCNP data files
    std::string clean_mcnp_command = "rm " + this->_run_directory + input_file_name + " " + this->_run_directory + output_file_name + " " + this->_run_directory + runtpe_name + " " + this->_run_directory + srctpe_name;
    exec(clean_mcnp_command);
    
    //create a symbolic link to the Doppler broadened cross sections
    std::string symbolic_link_command = "cd " + this->_run_directory + "; ln -s ../../../doppler-broadened-cs/otf*txt .";
    exec(symbolic_link_command);
    
    
    this->createMCNPOutputFile(run_title, input_file_name, delayed_neutrons);
    std::string command_line_log_file = "mcnp_run_log.txt";
    
    //Run the file
    #ifdef LAPTOP
    
    //We are just running MPI here
    std::string mcnp_path = "/media/chris/DoubleSpace/MCNP/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    std::string command = "cd " + this->_run_directory + "; mpirun -np  7 " + mcnp_path + " i=\"" + input_file_name + "\" o=\"" + output_file_name + "\" runtpe=\"" + runtpe_name + "\" srctp=\"" + srctpe_name + "\" | tee \"" + command_line_log_file + "\"";
    exec(command);
    
    #elif PRACTICE_CLUSTER
    
    std::string mcnp_path = "/share/apps/mcnp/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    
    //Delete the old log output by overwriting it with a blank string (we are using the log to check to see if the qsub is done)
    std::string remove_command = "cd " + this->_run_directory + "; echo \"\" > " + command_line_log_file;
    exec(remove_command);

    //Run Submission Script with the created MCNP file
    std::string qsub_command = "cd " + this->_run_directory + ";qsub -N " + this->_reactor->_run_name + " -pe orte 32 ../../../composite-fuel-submission-script.sh " + file_root + " " + command_line_log_file;
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

std::string ReactorMonteCarlo::getSingleCellCard(const Materials &material, const int &current_zone, int &cell_number )
{
    std::stringstream cell_card; 
    
     
    //MCNP reads temperature in MeV so we need a conversion factor for our programs K
    const Real MeVperK = 8.617e-11;
    //For now we will assume that the density stays constant regardless of temperature change to perserve conservation of mass
    Real density_derived_temperature = 400;
    //First is the density, second is the density derivative which we don't need. Divide by 1000 to convert from kg/m^3 to g/cm^3
    Real density = _reactor->_micro_sphere_geometry->_material_library.getDensityPair(material,density_derived_temperature,0).first/1000;

    for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
    {
        
        //Gather the temperature for the zone in this cell
        Real temperature;
        Real cell_volume;
        
        if( _cells_per_zone == 1) 
        {
            _reactor->_thermal_solver->getAverageTemperature(current_zone - 1, temperature, cell_volume);
        }
        else
        {
            _reactor->_thermal_solver->getCellTemperature(current_zone - 1, _cells_per_zone, current_cell_in_zone, temperature, cell_volume );
        }
    
        cell_volume *= (100 * 100 * 100);
        
        //we want to put our reflecting boundary condition here
        if( cell_number == 1 )
        {
            cell_card << " " << cell_number << " " << current_zone << " -" << density << " -"  << cell_number << " imp:n=1 TMP=" << temperature*MeVperK  << " $ " << getMaterialName(material) << " T = " << temperature << " K" << " Volume = " << cell_volume << " cm^3 " <<std::endl;
        }
        else
        {
            cell_card << " " << cell_number << " " << current_zone << " -" << density << " " << ( cell_number -1 )  << " " << " -"  << cell_number << " imp:n=1 TMP=" << temperature*MeVperK <<" $ " << getMaterialName(material) << " T = " << temperature << " K" << " Volume = " << cell_volume << " cm^3 " << std::endl;
        }

        cell_number++;
    }
    return cell_card.str();
}


std::string ReactorMonteCarlo::getCellCards()
{
    std::stringstream cell_cards;
    
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    int number_zones = geometry_data.size();   
    int current_cell_number = 1;
    
    //For each zone create a cell
    for( int index = 0; index < number_zones  ; index++ )
    {
        Materials material = geometry_data[index].first;
        int current_zone = index + 1;        
        //Note that the current_cell_number is passed by reference and its value is changed
        cell_cards << this->getSingleCellCard(material,current_zone, current_cell_number);                    
    }
    
    return cell_cards.str();
}

std::string ReactorMonteCarlo::getSurfaceCards()
{
    std::stringstream surface_cards;
    
    int surface_card_number = 1;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    int number_zones = geometry_data.size();   
            
    
    for( int index = 0; index < number_zones  ; index++ )
    {
        int current_zone = index + 1;
        
        if(this->_cells_per_zone == 1)
        {    
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
        else
        {
            Real last_zone_shell_radius = 0;
            
            if( index > 0)
            {
                last_zone_shell_radius = geometry_data[index -1 ].second * 100;
            }
            
            Real zone_shell_radius = geometry_data[index].second * 100; 
            
            for( int current_surface = 1; current_surface <= _cells_per_zone; current_surface++ )
            {
                //we want to put our reflecting boundary condition here
               if( current_zone == number_zones && current_surface == _cells_per_zone)
               {
                   surface_cards << "*";
               }
               else
               {
                   surface_cards << " ";
               }
               //cm sphere radius
               
               Real surface_radius = last_zone_shell_radius + (zone_shell_radius - last_zone_shell_radius ) * (static_cast<Real>(current_surface)/static_cast<Real>(_cells_per_zone));
               
               surface_cards << surface_card_number << " SPH 0 0 0 " << surface_radius << std::endl;
               surface_card_number++;
            }
            
            
        }
        
    }
    
    return surface_cards.str();
}


void ReactorMonteCarlo::createMCNPOutputFile(const std::string &run_title, const std::string &file_name,const bool &delayed_neutrons)
{
    
    std::stringstream mcnp_file;
   
    std::string cell_cards = this->getCellCards();
    std::string surface_cards = this->getSurfaceCards();
    std::string material_cards = this->getMaterialCards();
    int number_of_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles", 33);
    
    mcnp_file << run_title << std::endl;
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
    
    if(!delayed_neutrons)
    {
        mcnp_file << " TOTNU NO" << std::endl;
    }
    
    
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
