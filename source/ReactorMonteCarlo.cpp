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

ReactorMonteCarlo::ReactorMonteCarlo(InfiniteCompositeReactor* reactor,const Real &starting_k_eff, const  std::string &run_directory)
{
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _reactor = reactor;
    
    //Cells per zone splits each zone up into multiple cells
    _cells_per_zone = this->_reactor->_input_file_reader->getInputFileParameter("Cells Per Zone", 1 );  
    _number_zones = _reactor->_micro_sphere_geometry->_geometry.size();
    //Cells per zone splits each zone up into multiple cells
    _number_cpus = this->_reactor->_input_file_reader->getInputFileParameter("Number CPUs", 32 );  
    _k_eff_number_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles", 60 );  
    _beta_eff_number_cycles = this->_reactor->_input_file_reader->getInputFileParameter("Number of MCNP Cycles for Beta", 600 );  
    _calulate_beta_interval = this->_reactor->_input_file_reader->getInputFileParameter("Keff Calculation Per Beta Eff Calculation",0);
    _tally_cells = this->_reactor->_input_file_reader->getInputFileParameter("Tally Cells", false);
    _tally_energy_bins = this->_reactor->_input_file_reader->getInputFileParameter("Tally Energy Bins", 40);   
    _starting_k_eff = starting_k_eff;
    
    _current_beta_eff = -1;
    _current_beta_eff_sigma = -1;
    _virtual_k_eff_multiplier = -1;
    _number_of_keff_calculations = 0;
    
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const BetaSimulationResults &beta_results)
{
    _current_beta_eff = beta_results._beta;
    _current_beta_eff_sigma = beta_results._beta_sigma;
    this->updateCurrentValuesFromResults( beta_results._with_delayed_neutrons);
}

void ReactorMonteCarlo::updateCurrentValuesFromResults(const SimulationResults &results)
{
    if( _virtual_k_eff_multiplier < 0)
    {    
        _virtual_k_eff_multiplier = _starting_k_eff / results._k_eff;  
        _number_of_keff_calculations = 1;
    }
    
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
    std::string clean_mcnp_command = "rm " + this->_run_directory + input_file_name + " " + this->_run_directory + output_file_name + " " + this->_run_directory + srctpe_name + " " + this->_run_directory + mctal_name;
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
    
    #elif ACORITE_CLUSTER
    
    std::string submission_script = "acorite-submission-script.sh";
    
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
    }while(is_done == "");
    
    
    
    #endif
    
    //If we have the tally set take the results of the tally
    if(_tally_cells)
    {  
        TallyGroup* tally_group = TallyGroup::MCNPTallyGroupFactory( _run_directory + mctal_name, _number_zones, _cells_per_zone, _reactor->_transient_time);
        _tally_groups.push_back(tally_group);
        this->outputTalliesToFile(tally_group, _reactor->_results_directory + "tallydata.csv");
    }
    //Remove symbolic links to the Doppler broadened cross sections and the source tape
    //These are the largest files and we don't need them so remove them now
    std::string rm_symbolic_link_command = "cd " + this->_run_directory + "; rm otf*txt " + runtpe_name;
    exec(rm_symbolic_link_command);
    
    
    //Read the output file
    SimulationResults results = SimulationResults(output_file_name, this->_run_directory );
    return results;
        
}

std::vector< std::vector<Real> > ReactorMonteCarlo::getZoneCellRelativePowerDensity()
{
    std::vector<std::vector<Real>> cell_zone_power_densities = std::vector<std::vector<Real>>();
    
    for(int zone = 1; zone <= _number_zones; zone++)
    {
        std::vector<Real> cell_power_densities = std::vector<Real>();

        for( int cell = 1; cell <= _cells_per_zone; ++cell)
        {   
            if(! _tally_cells || _tally_groups.size() == 0)
            {
                
                if( cell == 1 && zone == 1)
                {
                    cell_power_densities.push_back(static_cast<Real>(1.0));
                }
                else
                {
                    cell_power_densities.push_back(static_cast<Real>(0.0));
                }
            }
            else
            {
                TallyGroup* latest_tally_group = this->_tally_groups.back();    
                Real power_density = latest_tally_group->_fission_tallies[ (zone-1)*_cells_per_zone + (cell-1) ]->_value;
                cell_power_densities.push_back(power_density);
            }
            
        }
        
        cell_zone_power_densities.push_back(cell_power_densities);
    }
    
    return cell_zone_power_densities;
}

std::string ReactorMonteCarlo::getMaterialCards()
{
    std::stringstream material_cards, otfdb_card;
        
    Real enrichment_fraction = _reactor->_input_file_reader->getInputFileParameter("Uranium Enrichment Fraction", static_cast<Real>(0.2) );

    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
     
       
    for( size_t index = 0; index < _number_zones  ; index++ )
    {
        size_t current_zone = index + 1;
    
        Materials material = geometry_data[index].first; 
        std::string material_card_entry;
        std::string doppler_card_entry;
        MaterialLibrary::getMcnpMaterialCard(material,current_zone,material_card_entry, doppler_card_entry, enrichment_fraction);
        
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
    Real density = MaterialLibrary::getDensityPair(material,density_derived_temperature,0).first/1000;

    for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
    {
        
        //Gather the temperature for the zone in this cell
        Real temperature;
        Real cell_volume;
        
        _reactor->_thermal_solver->getCellTemperature(current_zone - 1, _cells_per_zone, current_cell_in_zone, temperature, cell_volume );
           
        cell_volume *= (100 * 100 * 100);
        
        //we want to put our reflecting boundary condition here
        if( cell_number == 1 )
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << current_zone << " " << std::setw(9) << -density << " " << std::setw(3) <<                       "  " << std::setw(4)  << -cell_number << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) << getMaterialName(material) << " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
        }
        else
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << current_zone << " " << std::setw(9) << -density << " " << std::setw(3) << ( cell_number -1 )  << " " << std::setw(4)  << -cell_number << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) << getMaterialName(material) << " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
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
    
    tally_cards << " E0 0.000000001 " + std::to_string(_tally_energy_bins - 2) + "ILOG 10" << std::endl;
    tally_cards << " PRDMP   j j 1 		$write mctal file" << std::endl;
    
    return tally_cards.str(); 
}

std::string ReactorMonteCarlo::getSurfaceCards()
{
    std::stringstream surface_cards;
    
    int surface_card_number = 1;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    
    //std::string shape_type;        
    
    for( int index = 0; index < _number_zones  ; index++ )
    {
        int current_zone = index + 1;       
        Real last_zone_shell_radius = 0;

        if( index > 0)
        {
            last_zone_shell_radius = geometry_data[index -1 ].second * 100;
        }

        Real zone_shell_radius = geometry_data[index].second * 100; 

        for( int current_surface = 1; current_surface <= _cells_per_zone; current_surface++ )
        {
           Real surface_radius = last_zone_shell_radius + (zone_shell_radius - last_zone_shell_radius ) * (static_cast<Real>(current_surface)/static_cast<Real>(_cells_per_zone));
            
            //we want to put our reflecting boundary condition here
           if( current_zone == _number_zones && current_surface == _cells_per_zone)
           {
               Real box_edge = surface_radius*2;
               Real half_box_edge = surface_radius;
               surface_cards << "*" << surface_card_number << " BOX ";
               surface_cards << -half_box_edge << " " << -half_box_edge << " " << -half_box_edge << "  ";
               surface_cards << box_edge << " 0 0  ";
               surface_cards << "0 " << box_edge << " 0  ";
               surface_cards << "0 0  " << box_edge << std::endl;
           }
           else
           {
               surface_cards << " " << surface_card_number << " SPH 0 0 0 " << surface_radius << std::endl;
           }
           //cm sphere radius

           
           
           surface_card_number++;         
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
    
    mcnp_file << " KCODE 2000 1.5 3 " << number_of_cycles << "  $need at least 30 active cycles to print results" << std::endl;
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
    
    output_file << "Time [s], Zone, Cell, Value, Sigma";
    
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
