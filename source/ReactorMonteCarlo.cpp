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
#include <iomanip>
#include "InfiniteCompositeReactor.h"
#include "ReactorMonteCarlo.h"
#include "SimulationResults.h"
#include "TallyGroup.h"

ReactorMonteCarlo::~ReactorMonteCarlo() 
{
    for(std::size_t index = 0; index < _tally_groups.size(); ++index)
    {
        delete _tally_groups[index];
    }
}

ReactorMonteCarlo::ReactorMonteCarlo(Reactor* reactor, const  std::string &run_directory)
{
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _reactor = reactor;
    
    //Cells per zone splits each zone up into multiple cells
    _cells_per_zone = this->_reactor->_input_file_reader->getInputFileParameter("Cells Per Zone", static_cast<int>(1) );  
    _number_zones = _reactor->_micro_sphere_geometry->_geometry.size();
    //Cells per zone splits each zone up into multiple cells
    _number_cpus = this->_reactor->_input_file_reader->getInputFileParameter("Number CPUs", static_cast<int>(32) );  
    
    //Here we want to establish the number of cycles based on the 
    _k_eff_number_particles = this->_reactor->_input_file_reader->getInputFileParameter("Keff Number of Particles", static_cast<long>(33000) );  
    _beta_eff_number_particles = this->_reactor->_input_file_reader->getInputFileParameter("Beff Number of Particles", static_cast<long>(33000) );  
    _particles_per_cycle = this->_reactor->_input_file_reader->getInputFileParameter("Particles Per Cycle", static_cast<long>(1000) );
    
    _beta_eff_number_cycles = std::ceil(static_cast<Real>(_beta_eff_number_particles)/static_cast<Real>(_particles_per_cycle));
    _k_eff_number_cycles = std::ceil(static_cast<Real>(_k_eff_number_particles)/static_cast<Real>(_particles_per_cycle));
    
    if(_beta_eff_number_cycles < 33 || _k_eff_number_cycles < 33)
    {
        throw std::string("Need at least 33 cycles for MCNP to run correctly only " + std::to_string(_k_eff_number_cycles) + " cycles now available.");
    }
    
    _calulate_beta_interval = this->_reactor->_input_file_reader->getInputFileParameter("Keff Calculation Per Beta Eff Calculation",0);
    _tally_cells = this->_reactor->_input_file_reader->getInputFileParameter("Tally Cells", false);
    _tally_energy_bins = this->_reactor->_input_file_reader->getInputFileParameter("Tally Energy Bins", 40);     
    
    _current_beta_eff = -1;
    _current_beta_eff_sigma = -1;
    _number_of_keff_calculations = 0;
    _current_mc_exection_elapsed_time = -1;
    
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const BetaSimulationResults &beta_results)
{
    _current_beta_eff = beta_results._beta;
    _current_beta_eff_sigma = beta_results._beta_sigma;
    this->updateCurrentValuesFromResults( beta_results._with_delayed_neutrons);
    _current_mc_exection_elapsed_time = beta_results._elapsed_time;
    //Twice as many because beta calculation do do cycles
    _current_number_particles = _beta_eff_number_particles*2.0;
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const SimulationResults &results)
{
    _current_k_eff = results._k_eff;
    _current_prompt_neutron_lifetime = results._prompt_neutron_lifetime;
    _current_k_eff_sigma = results._k_eff_sigma;
    _current_prompt_neutron_lifetime_sigma = results._prompt_neutron_lifetime_sigma;
    _current_mc_exection_elapsed_time = results._elapsed_time;
    _current_number_particles = _k_eff_number_particles;
}


/**
 * Adjusts the k_effective to match the desired starting reactivity
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::updateAdjustedCriticalityParameters()
{
    if( _calulate_beta_interval <= 0 || _number_of_keff_calculations % _calulate_beta_interval == 0 || _current_beta_eff < 0 )
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
    SimulationResults k_eff_with_delayed = getRawCriticalityParameters( "k-effective-calc", _particles_per_cycle, _beta_eff_number_cycles, true);
    SimulationResults k_eff_without_delayed = getRawCriticalityParameters( "k-effective-no-delayed-calc", _particles_per_cycle, _beta_eff_number_cycles, false);    
    //Create the BetaSimulationResults Object
    BetaSimulationResults beta_results = BetaSimulationResults(k_eff_with_delayed,k_eff_without_delayed);
    return beta_results;
}


SimulationResults ReactorMonteCarlo::getRawKeff()
{
    //return the simulation results
    return this->getRawCriticalityParameters("k-effective-calc", this->_particles_per_cycle, _k_eff_number_cycles, true);    
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
SimulationResults ReactorMonteCarlo::getRawCriticalityParameters(const std::string &file_root, const int &particles_per_cycle, const int &number_cycles, const bool &delayed_neutrons)
{
    
    
    //create the MCNP input file
    std::string run_title = this->_reactor->_run_name + " Time = "  + std::to_string(this->_reactor->_thermal_solver->_current_time);
    std::string input_file_name = file_root + ".inp";
    std::string output_file_name =  file_root + ".out";
    std::string runtpe_name =  file_root + ".runtpe";
    std::string srctpe_name =  file_root + ".srctp";
    std::string mctal_name =  file_root + ".mctal";
    
    //clean up the previous MCNP data files
    std::string clean_mcnp_command = "rm " + this->_run_directory + input_file_name + " " + this->_run_directory + output_file_name + " " + this->_run_directory + srctpe_name + " " + this->_run_directory + mctal_name;
    exec(clean_mcnp_command);
    
    //create a symbolic link to the Doppler broadened cross sections
    std::string symbolic_link_otf_db_command = "cd " + this->_run_directory + "; ln -s ../../../doppler-broadened-cs/otf*txt .";
    exec(symbolic_link_otf_db_command);
    
    if(Reactor::_otf_sab)
    {
        //create a symbolic link to the Doppler broadened cross sections
        std::string symbolic_link_otf_sab_command = "cd " + this->_run_directory + "; ln -s ../../../otfsab/*00t .";
        exec(symbolic_link_otf_sab_command);
    }
    
    this->createMCNPOutputFile(run_title, input_file_name, particles_per_cycle, number_cycles, delayed_neutrons);
    std::string command_line_log_file = "mcnp_run_log.txt";
    
    std::time_t mc_start_time = std::time(nullptr);
    
    //Run the file
    #ifdef LAPTOP
    
    //We are just running MPI here
    std::string mcnp_path = "/media/chris/DoubleSpace/MCNP/MCNP_CODE/MCNP611/bin/mcnp6.mpi";
    std::string datapath_command = "export DATAPATH=/media/chris/DoubleSpace/MCNP/MCNP_DATA/";
    //std::string xsdir = "/media/chris/DoubleSpace/MCNP/MCNP_DATA/xsdir_mcnp6.1";
    std::string command = datapath_command + ";cd " + this->_run_directory + "; PATH=/home/chris/Programs/openmpi/bin:$PATH;LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/chris/Programs/openmpi/lib; mpirun -np  7 " + mcnp_path + " i=\"" + input_file_name + "\" o=\"" + output_file_name + "\" runtpe=\"" + runtpe_name + "\" srctp=\"" + srctpe_name + "\" mctal=\"" + mctal_name + "\" | tee \"" + command_line_log_file + "\"";
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
    
    #elif ACORITE_CLUSTER
    
    std::string submission_script = "acorite-submission-script.sh";
    
    #elif NEAMS

    std::string submission_script;
    
    if(Reactor::_otf_sab)
    {
        submission_script = "neams-pavlou-submission.sh";
    }
    else
    {
        submission_script = "neams-submission-script.sh";
    }
    
    #endif
    
    
    //Run Submission Script with the created MCNP file
    std::string qsub_command = "cd " + this->_run_directory + ";qsub -N " + this->_reactor->_run_name + "-t-" + std::to_string(static_cast<float>(_reactor->_transient_time) ) + " -pe orte " + std::to_string(_number_cpus) + " ../../../job-submission/" + submission_script + " " + file_root + " " + command_line_log_file;
    exec(qsub_command);
    //Constantly read the output file until it says mcrun done 
    std::string search_lock = "cd " + this->_run_directory + ";cat " + command_line_log_file + " | grep \"mcrun  is done\"";
    std::string is_done;
    
    do
    {   
        std::this_thread::sleep_for(std::chrono::milliseconds(5000));
        is_done = exec(search_lock);
    }
    while(is_done == "");
    
    
    #endif
    
    std::time_t mc_end_time = std::time(nullptr);    
    std::time_t mc_elapased_time = mc_end_time - mc_start_time;
    
    //If we have the tally set take the results of the tally
    if(_tally_cells)
    {  
        TallyGroup* tally_group = TallyGroup::MCNPTallyGroupFactory( _run_directory + mctal_name, _number_zones, _cells_per_zone, _reactor->_transient_time);
        _tally_groups.push_back(tally_group);
        this->outputTalliesToFile(tally_group, _reactor->_results_directory + "tallydata.csv");
    }
    //Remove symbolic links to the Doppler broadened cross sections and the source tape
    //These are the largest files and we don't need them so remove them now
    std::string rm_symbolic_link_command = "cd " + this->_run_directory + "; rm otf*txt " + runtpe_name + "; rm *.00t";
    exec(rm_symbolic_link_command);
    
    
    //Read the output file
    SimulationResults results = SimulationResults(output_file_name, this->_run_directory, mc_elapased_time );
    return results;
        
}

std::vector< std::vector<Real> > ReactorMonteCarlo::getZoneCellRelativePowerDensity()
{
    std::vector<std::vector<Real>> cell_zone_power_densities = std::vector<std::vector<Real>>();
    
    Real maximum_power_density = -1;
    
    //Go through each zone and each cell in each zone 
    for(int zone = 1; zone <= _number_zones; zone++)
    {
        std::vector<Real> cell_power_densities = std::vector<Real>();

        for( int cell = 1; cell <= _cells_per_zone; ++cell)
        {   
            Real power_density;
            
            //If no tallies set up so that the first zone is powered
            if(! _tally_cells || _tally_groups.size() == 0)
            {
                
                if( cell == 1 && zone == 1)
                {
                    power_density = 1.0;
                }
                else
                {
                    power_density = 0.0;
                }
            }
            //Otherwise use tallies
            else
            {
                TallyGroup* latest_tally_group = this->_tally_groups.back();    
                power_density = latest_tally_group->_fission_tallies[ (zone-1)*_cells_per_zone + (cell-1) ]->_value;               
            }
            
            if(power_density > maximum_power_density)
            {
                maximum_power_density = power_density;
            }
            
            cell_power_densities.push_back(power_density);            
        }
        
        cell_zone_power_densities.push_back(cell_power_densities);
    }
    
    for(int zone = 0; zone < _number_zones; ++zone)
    {
        for( int cell = 0; cell < _cells_per_zone; ++cell)
        {               
            cell_zone_power_densities[zone][cell] = cell_zone_power_densities[zone][cell]/maximum_power_density;
        }
     }
    
    return cell_zone_power_densities;
}

std::string ReactorMonteCarlo::getTallyCards()
{
    std::stringstream tally_cards;
    
    
    tally_cards << "c Tally Cards" << std::endl;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    int number_zones = this->_number_zones;   
    int cell_number = 1;
    //Create the Fission and Flux Tallies
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
    
    //Create the Absorption Tallies
    int tally_number = cell_number;
    //cell number couts from 1 to n_zones * n_cells
    cell_number = 1;
    
    //For each zone create a cell
    for( int current_zone = 1; current_zone <= number_zones  ; current_zone++ )
    {
        for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
        {
            //Here we are creating our absorption tally
            tally_cards << " F" + std::to_string(tally_number) << "4:n " << std::to_string(cell_number) << "       $ absorption rate tally" << std::endl;
            tally_cards << " FM" + std::to_string(tally_number) << "4 -1 " << current_zone << " -6:-2  $ multiplier to set the absorption rate tally counter" << std::endl;
            cell_number++;
            tally_number++;
        }
    }
    
    //cell number couts from 1 to n_zones * n_cells
    cell_number = 1;
    
    //For each zone create a cell
    for( int current_zone = 1; current_zone <= number_zones  ; current_zone++ )
    {
        for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
        {
            //Here we are creating our fission number tally
            tally_cards << " F" + std::to_string(tally_number) << "4:n " << std::to_string(cell_number) << "       $ fission rate tally" << std::endl;
            tally_cards << " FM" + std::to_string(tally_number) << "4 -1 " << current_zone << " -2  $ multiplier to set the fission rate tally counter" << std::endl;
            cell_number++;
            tally_number++;
        }
    }
    
    //cell number couts from 1 to n_zones * n_cells
    cell_number = 1;
    
    for( int current_zone = 1; current_zone <= number_zones  ; current_zone++ )
    {
        for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
        {
            //Here we are creating our absorption tally
            tally_cards << " F" + std::to_string(tally_number) << "4:n " << std::to_string(cell_number) << "       $ capture rate tally" << std::endl;
            tally_cards << " FM" + std::to_string(tally_number) << "4 -1 " << current_zone << " -6  $ multiplier to set the capture rate tally counter" << std::endl;
            cell_number++;
            tally_number++;
        }
    }
    
    
    tally_cards << " E0 0.000000001 " + std::to_string(_tally_energy_bins - 2) + "ILOG 10" << std::endl;
    tally_cards << " PRDMP   j j 1 		$write mctal file" << std::endl;
    
    return tally_cards.str(); 
}

void ReactorMonteCarlo::createMCNPOutputFile(const std::string &run_title, const std::string &file_name,const int &paticles_per_cycle,const int &number_of_cycles, const bool &delayed_neutrons)
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
    
    if(Reactor::_otf_sab)
    {
        mcnp_file << " xs1 grph.00t 11.898000 grph.00t 0 1 1 1409891 0 0 0.000E-00" << std::endl;
    }
    
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
    
    mcnp_file << " KCODE " << paticles_per_cycle << " 1.5 3 " << number_of_cycles << "  $need at least 30 active cycles to print results" << std::endl;
    mcnp_file << " KSRC 0 0 0" << std::endl;
    mcnp_file << " print" << std::endl;
    mcnp_file << "c end data" << std::endl;
    
    std::ofstream mcnp_input_file;
    mcnp_input_file.open( this->_run_directory + file_name,std::ios::out);
    mcnp_input_file << mcnp_file.str();
    mcnp_input_file.close();
    
    std::cout<<mcnp_file.str();
}

void ReactorMonteCarlo::createTallyOutputFile(std::string tally_file_name)
{
    std::ofstream output_file;
    output_file.open( tally_file_name, std::ios::out);
    
    output_file << "Time [s],Name,Zone,Cell,Value,Sigma";
    
    for(size_t index = 0; index < _tally_energy_bins; index++ )
    {
        output_file << ",Energy-" << index <<" [MeV],value-" << index << ",sigma-" << index;
    }
    
    output_file << std::endl;    
    output_file.close();
}

void ReactorMonteCarlo::outputTalliesToFile(TallyGroup* tally_group, std::string file_path)
{
    if(!file_exists(file_path))
    {
        this->createTallyOutputFile(file_path);
    }
    
    tally_group->outputToFile(file_path);
}

Real ReactorMonteCarlo::getPowerPeaking()
{
    if(this->_tally_cells)
    {
        std::vector<std::vector<Real>> cell_zone_power_densities = this->getZoneCellRelativePowerDensity();
        Real min_value = 1;
        
        for(int zone = 0; zone < _number_zones; ++zone)
        {
            for( int cell = 0; cell < _cells_per_zone; ++cell)
            {               
                if( cell_zone_power_densities[zone][cell] < min_value  &&  cell_zone_power_densities[zone][cell] > 0)
                {
                    min_value = cell_zone_power_densities[zone][cell];
                }
            }
         }
        
        return 1.0/min_value;
    }
    else
    {
        return -1;
    }
}