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


#include "MicroGeometry.h"
#include <sstream>
#include <math.h>
#include "InputDataFunctions.h"
#include <chrono>
#include <thread>
#include "InfiniteCompositeReactor.h"
#include "ReactorMonteCarlo.h"
#include "SimulationResults.h"

ReactorMonteCarlo::ReactorMonteCarlo() {}

ReactorMonteCarlo::ReactorMonteCarlo(InfiniteCompositeReactor* reactor,const Real &starting_k_eff, const  std::string &run_directory)
{
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _reactor = reactor;
    
    //Cells per zone splits each zone up into multiple cells
    _cells_per_zone = this->_reactor->_input_file_reader->getInputFileParameter("Cells Per Zone", 1 );  
    
    //Cells per zone splits each zone up into multiple cells
    _number_cpus = this->_reactor->_input_file_reader->getInputFileParameter("Number CPUs", 32 );  
    _k_eff_number_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles", 60 );  
    _beta_eff_number_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles for Beta", 600 );  
    _calulate_beta_interval = this->_reactor->_input_file_reader->getInputFileParameter("Keff Calculation Per Beta Eff Calculation",0);
    _tally_cells = this->_reactor->_input_file_reader->getInputFileParameter("Tally Cells", false);
   
    
    BetaSimulationResults beta_results = getRawKeffAndBetaEff();
    
    
    _virtual_k_eff_multiplier = starting_k_eff / beta_results._with_delayed_neutrons._k_eff;  
    this->updateCurrentValuesFromResults(beta_results);
    this->_number_of_keff_calculations = 1;
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const BetaSimulationResults &beta_results)
{
    _current_beta_eff = beta_results._beta;
    _current_beta_eff_sigma = beta_results._beta_sigma;
    this->updateCurrentValuesFromResults( beta_results._with_delayed_neutrons);
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const SimulationResults &results)
{
    _current_k_eff = results._k_eff * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime = results._prompt_neutron_lifetime;
    _current_k_eff_sigma = results._k_eff_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = results._prompt_neutron_lifetime_sigma;
}


/**
 * Adjusts the k_effective to match the desired starting reactivity
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::updateAdjustedCriticalityParameters()
{
    if( _calulate_beta_interval <= 0 || _number_of_keff_calculations % _calulate_beta_interval == 0 )
    {
        BetaSimulationResults beta_results = this->getRawKeffAndBetaEff();
        this->updateCurrentValuesFromResults(beta_results);
    }
    else
    {
        SimulationResults results = this->getRawKeff();
        this->updateCurrentValuesFromResults(results);
    }
    
    _number_of_keff_calculations++;
}

BetaSimulationResults ReactorMonteCarlo::getRawKeffAndBetaEff()
{
    //Run the simulation to get the results for both delated and non delayed to do the beta calculation
    SimulationResults k_eff_with_delayed = getRawCriticalityParameters( "keffective-calc", _beta_eff_number_cycles, true);
    SimulationResults k_eff_without_delayed = getRawCriticalityParameters( "keffective-no-delayed-calc", _beta_eff_number_cycles, false);    
    //Create the BetaSimulationResults Object
    BetaSimulationResults beta_results = BetaSimulationResults(k_eff_with_delayed,k_eff_without_delayed);
    return beta_results;
}


SimulationResults ReactorMonteCarlo::getRawKeff()
{
    //return the simulation results
    return this->getRawCriticalityParameters("k-effective-calc", _k_eff_number_cycles, true);    
}

/**
 * Create the input file, run it, and grab the k-effective from the MCNP output file
 * 
 * @param file_root
 * @param k_eff
 * @param k_eff_sigma
 * @param prompt_removal_lifetime
 * @param prompt_removal_lifetime_sigma
 * @param number_cycles
 * @param delayed_neutrons
 */
SimulationResults ReactorMonteCarlo::getRawCriticalityParameters(const std::string &file_root, const int &number_cycles, const bool &delayed_neutrons)
{
    //create the MCNP input file
    std::string run_title = this->_reactor->_run_name + " Time = "  + std::to_string(this->_reactor->_thermal_solver->_current_time);
    std::string input_file_name = file_root + ".inp";
    std::string output_file_name =  file_root + ".out";
    std::string runtpe_name =  file_root + ".runtpe";
    std::string srctpe_name =  file_root + ".srctp";
    std::string mctal_name =  file_root + ".mctal";
    
    //clean up the previous MCNP data files
    std::string clean_mcnp_command = "rm " + this->_run_directory + input_file_name + " " + this->_run_directory + output_file_name + " " + this->_run_directory + runtpe_name + " " + this->_run_directory + srctpe_name + " " + this->_run_directory + mctal_name;
    exec(clean_mcnp_command);
    
    //create a symbolic link to the Doppler broadened cross sections
    std::string symbolic_link_command = "cd " + this->_run_directory + "; ln -s ../../../doppler-broadened-cs/otf*txt .";
    exec(symbolic_link_command);
    
    
    this->createMCNPOutputFile(run_title, input_file_name, number_cycles, delayed_neutrons);
    std::string command_line_log_file = "mcnp_run_log.txt";
    
    //Run the file
    #ifdef LAPTOP
    
    //We are just running MPI here
    std::string mcnp_path = "/media/chris/DoubleSpace/MCNP/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    std::string command = "cd " + this->_run_directory + "; PATH=/home/chris/Programs/openmpi/bin:$PATH;LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/chris/Programs/openmpi/lib; mpirun -np  7 " + mcnp_path + " i=\"" + input_file_name + "\" o=\"" + output_file_name + "\" runtpe=\"" + runtpe_name + "\" srctp=\"" + srctpe_name + "\" mctal=\"" + mctal_name + "\" | tee \"" + command_line_log_file + "\"";
    exec(command);
    
    #else
    
    //Delete the old log output by overwriting it with a blank string (we are using the log to check to see if the qsub is done)
    std::string remove_command = "cd " + this->_run_directory + "; echo \"\" > " + command_line_log_file;
    exec(remove_command);

    #ifdef PRACTICE_CLUSTER
    
    std::string submission_script = "practice-cluster-submission-script.sh";
    
    #elif ANTAL
    
    std::string submission_script = "antal-submission-script.sh";
    
    #elif MC_CLUSTER
    
    std::string submission_script = "mc-submission-script.sh";
    
    #endif
    
    
    //Run Submission Script with the created MCNP file
    std::string qsub_command = "cd " + this->_run_directory + ";qsub -N " + this->_reactor->_run_name + " -pe orte " + std::to_string(_number_cpus) + " ../../../job-submission/" + submission_script + " " + file_root + " " + command_line_log_file;
    exec(qsub_command);
    //Constantly read the output file until it says mcrun done 
    std::string search_lock = "cd " + this->_run_directory + ";cat " + command_line_log_file + " | grep \"mcrun  is done\"";
    std::string is_done;
    
    do
    {   
        std::this_thread::sleep_for(std::chrono::milliseconds(5000));
        is_done = exec(search_lock);
    }while(is_done == "");
    
    
    
    #endif
    

    //Read the output file
    SimulationResults results = SimulationResults(output_file_name, this->_run_directory );
    return results;
        
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
    
    std::string U238_cs = "92238.80c";
    std::string U235_cs = "92235.80c";    

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

std::string ReactorMonteCarlo::getTallyCards()
{
    std::stringstream tally_cards;
    
    
    tally_cards << "c Tally Cards" << std::endl;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    int number_zones = geometry_data.size();   
    int cell_number = 1;
    
    //For each zone create a cell
    for( int current_zone = 1; current_zone <= number_zones  ; current_zone++ )
    {
        for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
        {
            tally_cards << " F" + std::to_string(cell_number) << "7:n " << std::to_string(cell_number) << "   $cell tally" << std::endl;
            tally_cards << " F" + std::to_string(cell_number) << "4:n " << std::to_string(cell_number) << "   $fission energy deposition tally" << std::endl;
            cell_number++;
        }
    }
    
    tally_cards << " E0 0.00000001 45ILOG 10" << std::endl;
    tally_cards << " PRDMP   j j 1 		$write mctal file" << std::endl;
    
    return tally_cards.str(); 
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


void ReactorMonteCarlo::createMCNPOutputFile(const std::string &run_title, const std::string &file_name,const int &number_of_cycles, const bool &delayed_neutrons)
{
    
    std::stringstream mcnp_file;
   
    std::string cell_cards = this->getCellCards();
    std::string surface_cards = this->getSurfaceCards();
    std::string material_cards = this->getMaterialCards();
    
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
    
    if(_tally_cells)
    {
        std::string tally_cards = this->getTallyCards();
        mcnp_file << tally_cards;
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
