/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ResultsStorage.cpp
 * Author: chris
 * 
 * Created on December 30, 2015, 3:01 PM
 */

#include <string>
#include <tuple>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include <cmath>
#include <iomanip>
#include <ctime>
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InputFileParser.h"

const std::time_t InfiniteCompositeReactor::_simulation_start_time = std::time(nullptr);;

InfiniteCompositeReactor::InfiniteCompositeReactor(const std::string &input_file_name ) 
{   
    
    //Check to make sure the old dire
    if( ! file_exists(input_file_name) )
    {
        std::cerr << "Input File: " + input_file_name + " doesn't exist";
        throw 1;

    }
    
    //Initialize the input file reader
    this->_input_file_reader = new InputFileParser( input_file_name );
    
    //Create the Results Folder
    time_t run_identification_number = std::time(nullptr);
    
    std::string default_run_name = "Unnamed-Run";
    _run_name = _input_file_reader->getInputFileParameter("Run Name", default_run_name);
    
    _results_directory =  "results/" + _run_name + "-" + std::to_string(run_identification_number) + "/";    
    
    _data_file = "datafile.csv";
    std::string folder_command = "mkdir -p " + _results_directory;
    exec( folder_command );
    
    std::string copy_input_file_command = "cp " + input_file_name + " " + _results_directory + "input_file.inp";
    exec( copy_input_file_command );
    
    _monte_carlo_number_iterations = 0;
    _transient_time = -1;
    
    initializeInifiniteCompositeReactorProblem();
}

/**
 * This constructor loads the conditions of a previously run program and continues it until the new end time
 * @param old_results_folder  the folder with the old files in it
 * @param new_end_time        the new end time (must be longer than the old one
 */
InfiniteCompositeReactor::InfiniteCompositeReactor(const std::string old_results_folder, Real new_end_time)
{
    
    //Check to make sure the old results folder exists 
    if( ! file_exists(old_results_folder) )
    {
        std::cerr << "Results Directory: " + old_results_folder + " Doesn't exist";
        throw 1;

    }

    //Then check for the files that are needed
    std::string old_input_file = old_results_folder + "/input_file.inp";
    std::string old_data_file = old_results_folder + "/datafile.csv";
    std::string old_temeprature_file = old_results_folder + "/temperature-data.csv";

    if( ! file_exists(old_input_file) || ! file_exists(old_data_file) || ! file_exists(old_temeprature_file))
    {
        std::cerr << "Missing Data Files " + old_input_file + " " + old_data_file + old_temeprature_file;
        throw 1;
    }
    
    //Initialize the input file reader
    this->_input_file_reader = new InputFileParser( old_input_file );
    
    //Create the Results Folder
    time_t run_identification_number = std::time(nullptr);
    
    std::string default_run_name = "Unnamed-Run";
    _run_name = _input_file_reader->getInputFileParameter("Run Name", default_run_name);
    
    _results_directory =  "results/" + _run_name + "-" + std::to_string(run_identification_number) + "/";    
    
    _data_file = "datafile.csv";
    std::string folder_command = "mkdir -p " + _results_directory;
    exec( folder_command );
    
    //Copy the files to the new directory
    std::string copy_input_file_command = "cp " + old_input_file + " " + _results_directory + "input_file.inp";
    exec( copy_input_file_command );
    
    std::string copy_old_data_file = "cp " + old_data_file + " " + _results_directory + "input_file.inp";
    exec( copy_old_data_file );
    
    std::string copy_temperature_file_command = "cp " + old_temeprature_file + " " + _results_directory + "input_file.inp";
    exec( copy_temperature_file_command );
    
    _monte_carlo_number_iterations = 0;
    _transient_time = 0;
    
    //initializeInifiniteCompositeReactorProblem();
}

void InfiniteCompositeReactor::simulate()
{
    this->_transient_time = 0;
    //Save the data for the initial timestep
    
    
    //Simulate the transient the outer loop is the monte carlo simulation
    for( ; _transient_time < _end_time; _transient_time += _inner_time_step)
    {
        this->monteCarloTimeStepSimulationProcessing();
        
        if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Time)
        {
            this->timeIterationInnerLoop();        
        }
        else if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Temperature )
        {
            this->temperatureIterationInnerLoop();
        }
        else
        {
            throw "Monte Carlo Iteration Type Unknown";
        }
    }
    
    //Save data and create the graphs for the post simulation data
    this->postSimulationProcessing();
}

void InfiniteCompositeReactor::monteCarloTimeStepSimulationProcessing()
{
    //Set the current power
    Real current_power = this->_kinetics_model->_current_power;

    //Gather the parameters from the monte carlo model 
    //The Monte Carlo model is run on the outer loop
    Real prompt_removal_lifetime = _monte_carlo_model->_current_prompt_neutron_lifetime;
    Real k_eff = _monte_carlo_model->_current_k_eff; 
    Real k_eff_sigma = _monte_carlo_model->_current_k_eff_sigma;
    Real beta_eff = _monte_carlo_model->_current_beta_eff;
    Real beta_eff_sigma = _monte_carlo_model->_current_beta_eff_sigma;        
    Real lambda = prompt_removal_lifetime/k_eff;
    Real lambda_sigma = _monte_carlo_model->_current_prompt_neutron_lifetime_sigma/k_eff;

    //Create the data structures for each time step's storage
    std::tuple<Real,Real,Real> k_effective_data = std::make_tuple( _transient_time , k_eff, k_eff_sigma );
    _k_eff_record.push_back( k_effective_data );

    std::tuple<Real,Real,Real> beta_effective_data = std::make_tuple( _transient_time , beta_eff, beta_eff_sigma );
    _beta_eff_record.push_back( beta_effective_data );

    std::tuple<Real,Real,Real> prompt_removal_lifetime_pair = std::make_tuple(_transient_time, lambda, lambda_sigma);
    _prompt_life_time_record.push_back(prompt_removal_lifetime_pair);

    Real reactivity = (k_eff - 1.0)/k_eff;
    Real reactivity_uncertainty = std::abs( reactivity ) * ( std::abs(k_eff_sigma / ( k_eff - 1)) + std::abs( k_eff_sigma / k_eff ) );
    Real reactivity_pcm = 10000.0 * reactivity;
    Real reactivity_pcm_uncertainty = 10000.0 * reactivity_uncertainty; 
    Real reactivity_cents = reactivity / beta_eff;
    Real reactivity_cents_uncertainty = std::abs( reactivity_cents ) * ( std::abs( beta_eff_sigma / beta_eff) + std::abs( reactivity_uncertainty / reactivity ));

    std::tuple<Real,Real,Real> reactivity_pcm_pair = std::make_tuple(_transient_time, reactivity_pcm, reactivity_pcm_uncertainty);
    _reactivity_pcm_record.push_back(reactivity_pcm_pair);        

    std::tuple<Real,Real,Real> reactivity_cents_pair = std::make_tuple(_transient_time, reactivity_cents, reactivity_cents_uncertainty);
   _reactivity_cents_record.push_back(reactivity_cents_pair);

    Real hottest_temperature = this->_thermal_solver->_solution[0]; 
    Real gamma = this->_monte_carlo_model->_virtual_k_eff_multiplier;
    
    this->saveCurrentData(_transient_time, current_power, k_eff, k_eff_sigma, lambda, lambda_sigma, beta_eff, beta_eff_sigma, hottest_temperature, gamma);
    
    MicroSolution solution = this->_thermal_solver->getCurrentMicrosolution();
    
    solution.plot( this->_results_directory + "solution-" + std::to_string( _transient_time )  + ".png", 800, 2000);
    _plot_solutions.push_back(solution); 
    
    //Save the current thermal solution to the output file
    std::vector<MicroSolution> current_solution = { solution };
    MicroSolution::saveSolutions( current_solution, this->_results_directory );
    
    this->_thermal_solver->logPowerDensity();
        
}

void InfiniteCompositeReactor::postSimulationProcessing()
{
    //Post Processing graph creation
    //MicroSolution::saveSolutions( _plot_solutions, this->_results_directory );
    MicroSolution::plotSolutions( _plot_solutions, 6 , this->_results_directory + "solutions-graph.png");
    PythonPlot::plotData(      _power_record,            "Time [s]", "Power Density [W/m^3]",      "", "Power vs. Time",                   this->_results_directory + "power-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _prompt_life_time_record, "Time [s]", "Prompt Neutron Lifetime [s]","", "Prompt Neutron Lifetime vs. Time", this->_results_directory + "prompt-neutron-lifetime-graph.png", {0, _end_time} );
    PythonErrorPlot::plotData( _k_eff_record,            "Time [s]", "K effective",    "",             "K-eff vs. Time",                   this->_results_directory + "k-eff-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_cents_record, "Time [s]", "Reactivity [Dollars]",    "",    "Reactivity vs. Time",              this->_results_directory + "reactivity-cents-graph.png",        {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_pcm_record,   "Time [s]", "Reactivity [pcm]",    "",        "Reactivity vs. Time",              this->_results_directory + "reactivity-pcm-graph.png",          {0, _end_time} );
    PythonErrorPlot::plotData( _beta_eff_record,         "Time [s]", "Beta effective",    "",          "Beta-eff vs. Time",                this->_results_directory + "beta-eff-graph.png",                {0, _end_time} );
    PythonPlot::plotData(      _delayed_record,          "Time [s]", "Delayed Precursors",         {}, "Keff vs. Delayed Precursors",      this->_results_directory + "delayed-precursors.png",            {0, _end_time} );
    PythonPlot::createPlots();        
}

bool InfiniteCompositeReactor::significantTemperatureDifference(MicroSolution* comparison)
{
     MicroSolution* reference = &_plot_solutions.back();
     
     MicroSolution differences = MicroSolution::temperatureDifference( comparison, reference);
     
     Real max_difference = 0;
     Real average_difference = 0;
     
     for(std::size_t ii=0; ii < differences.size(); ++ii)
     {
         if(max_difference < differences._solution[ii])
         {
             max_difference = differences._solution[ii];
         }
         
         if(ii !=0)
         {
            average_difference += differences._solution[ii] * ( std::pow(differences._grid[ii],3) - std::pow(differences._grid[ii - 1],3) );         
         }
         else
         {
             average_difference += differences._solution[ii] * ( std::pow(differences._grid[ii],3) );         
         }
     }     
     
     average_difference /= std::pow(differences._grid.back(),3);
     
     if( max_difference > _maximum_allowed_temperature_difference || average_difference > _maximum_allowed_temperature_difference/3 )
     {
         return true;
     }
     else
     {
         return false;
     }     
}

void InfiniteCompositeReactor::temperatureIterationInnerLoop()
{
    Real last_reported_time = 0;
    
    std::vector<std::vector<Real>> tally_power_density_map = _monte_carlo_model->getZoneCellRelativePowerDensity();
        
    //This loop iterates the kinetics and thermal model data
    for( _inner_time_step = 0 ; true /*infinite loop this and exit with a break*/ ; _inner_time_step += _kinetics_thermal_sync_time_step)
    {
        //Solve the kinetics model
        Real current_power = _kinetics_model->solveForPower(_kinetics_thermal_sync_time_step);
        std::vector<Real> power_distribition;
        
       
        if(_monte_carlo_model->_tally_cells)
        {
            //Get a vector representation of the radial power distribution
            power_distribition = _thermal_solver->getTallyBasedRepresentativeKernelPowerDistribution(tally_power_density_map, current_power );
        }
        else
        {
            power_distribition = _thermal_solver->getRespresentativePowerDistribution( current_power );
        }
        //Get the thermal solution
        MicroSolution current_solution = _thermal_solver->solve( _kinetics_thermal_sync_time_step, power_distribition); 

        
        //Record the power and delayed records on an intermediate time scale
        if( ( _inner_time_step - last_reported_time )  >= _power_and_delayed_neutron_record_time_step )
        {
            last_reported_time = _inner_time_step;

            Real absolute_time = _transient_time + _inner_time_step;

            //Everytime the MC is recalculated store this data
            std::pair<Real,Real> power_entry = { absolute_time, current_power };
            _power_record.push_back(power_entry);

            std::pair<Real,std::vector<Real>> delayed_entry = { absolute_time, _kinetics_model->_delayed_precursors };
            _delayed_record.push_back(delayed_entry);
        }
        
        if( this->significantTemperatureDifference(&current_solution))
        {
            break;
        }

    }

    
    _monte_carlo_model->updateAdjustedCriticalityParameters();
    _monte_carlo_time_iteration = _inner_time_step;   
    _monte_carlo_number_iterations++;
}

void InfiniteCompositeReactor::timeIterationInnerLoop()
{
    Real last_reported_time = 0;
        
    //This loop iterates the kinetics and thermal model data
    for( _inner_time_step = 0 ; _inner_time_step < _monte_carlo_time_iteration ; _inner_time_step += _kinetics_thermal_sync_time_step)
    {
        //Solve the kinetics model
        Real current_power = _kinetics_model->solveForPower(_kinetics_thermal_sync_time_step);

        //Get a vector representation of the radial power distribution
        std::vector<Real> power_distribition = _thermal_solver->getRespresentativePowerDistribution( current_power /* current shape */ );

        //Get the thermal solution
        _thermal_solver->solve( _kinetics_thermal_sync_time_step, power_distribition); 

        //Add the current time steps solution to the record (we can record all of them every 20 us let usually do every millisecond)            
        if( ( _inner_time_step - last_reported_time )  >= _power_and_delayed_neutron_record_time_step )
        {
            last_reported_time = _inner_time_step;

            Real absolute_time = _transient_time + _inner_time_step;

            //Everytime the MC is recalculated store this data
            std::pair<Real,Real> power_entry = { absolute_time, current_power };
            _power_record.push_back(power_entry);

            std::pair<Real,std::vector<Real>> delayed_entry = { absolute_time, _kinetics_model->_delayed_precursors };
            _delayed_record.push_back(delayed_entry);
        }

    }

    this->monteCarloTimeStepSimulationProcessing();
    
    //if there is still enough time left to do another monte carlo time iteration
    if( _transient_time + _inner_time_step < _end_time)
    {
        Real last_k_eff = _monte_carlo_model->_current_k_eff;

        
        _monte_carlo_model->updateAdjustedCriticalityParameters();

        Real current_k_eff =  _monte_carlo_model->_current_k_eff;

        Real difference = current_k_eff - last_k_eff;

        Real k_eff_change = std::abs( difference);


        Real k_eff_sigma = _monte_carlo_model->_current_k_eff_sigma;
        Real min_k_eff_change = k_eff_sigma*2/3;
        Real max_k_eff_change = k_eff_sigma*2.5;

        //We need smaller time steps
        if( max_k_eff_change < k_eff_change )
        {
            _monte_carlo_time_iteration *= 0.5;
        }
        //we can get away with bigger time steps
        else if( min_k_eff_change > k_eff_change )
        {
            _monte_carlo_time_iteration *= 1.5;
        }       
    }
    
     _monte_carlo_number_iterations++;
    
}

InfiniteCompositeReactor::~InfiniteCompositeReactor()
{
    delete _micro_sphere_geometry;
    delete _thermal_solver;
    delete _kinetics_model;
    delete _monte_carlo_model;
    delete _input_file_reader;
}
    
void InfiniteCompositeReactor::initializeInifiniteCompositeReactorProblem()
{
    
    std::string mc_recalc_type = _input_file_reader->getInputFileParameter(std::string("Monte Carlo Recalculation Type"), std::string("Temperature"));
    _monte_carlo_reclaculation_type = InfiniteCompositeReactor::getRecalculationType(mc_recalc_type);
    
    if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Temperature )
    {
        _maximum_allowed_temperature_difference = _input_file_reader->getInputFileParameter("MC Recalc Delta T", static_cast<Real>(20.0));
        _monte_carlo_time_iteration = -1;
    }
    else if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Time )
    {
        _monte_carlo_time_iteration =  _input_file_reader->getInputFileParameter("Monte Carlo Recalculation Timestep", static_cast<Real>(0.01) );  //How often to calculate keff and the prompt neutron lifetime
    }
    
    //This command establishes the logging of python plot commands. Useful for debugging
    PythonPlot::_log_file = this->_results_directory + "graph_log.log";
    
    //Default geometry sizes
    Real default_sphere_outer_radius = 2e-3;  //meters
    Real default_fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { default_fuel_kernel_outer_radius, default_sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    //Grab geometry from input file
    std::vector<Dimension> dimensions = _input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =   _input_file_reader->getInputFileParameter("Materials", default_materials);
    
    //Create a new geometry
    this->_micro_sphere_geometry = new MicroGeometry(materials, dimensions);    
          
    //Define the heat transfer settings
    Real initial_power_density =  _input_file_reader->getInputFileParameter("Starting Power Density",static_cast<Real>(200e6) ); // W/m^3 averaged over the entire micro sphere
    Real initial_outer_shell_temperature = 800;//_input_file.getInputFileParameter("Kernel Outer Temperature",800); // Kelvin
    
    //Set up a boundary condition 
    MicroCellBoundaryCondition* fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(initial_outer_shell_temperature);
    //Then crate the MicoCell thermal solver, set the boundary condition and then iterate steady state condition 
    this->_thermal_solver = new MicroCell(this, initial_outer_shell_temperature);
    this->_thermal_solver->setBoundaryCondition(fixed_temperature_boundary_condition);
    std::vector<Real> homogenous_power_density = this->_thermal_solver->getRespresentativePowerDistribution(initial_power_density);
    std::vector<MicroSolution> plot =  this->_thermal_solver->iterateInitialConditions(homogenous_power_density);
    MicroSolution::plotSolutions(plot,0, this->_results_directory + "initial-pre-tally-solve.png");
    
    //Define the Monte Carlo Parameters
    Real starting_k_eff =  _input_file_reader->getInputFileParameter("Starting K-eff",static_cast<Real>(1.01) );   
    this->_monte_carlo_model = new ReactorMonteCarlo(this, starting_k_eff, this->_results_directory + "run/");   
    
    if(_monte_carlo_model->_tally_cells)
    {
        this->solveForSteadyStatePowerDistribution(homogenous_power_density,initial_power_density);            
    }    
        
    this->_monte_carlo_model->updateAdjustedCriticalityParameters();
     
    
    //Reset the time to zero and then grab the current solution for the Problem
    _thermal_solver->_current_time = 0;    
     //Add the starting MicroSolution
    _plot_solutions.push_back( _thermal_solver->getCurrentMicrosolution() );
        
    //Get the steady state heat flux and set it as the boundary condition
    Real outer_boundary_heat_flux = sphere_volume(this->_thermal_solver->_mesh->_outer_radius[this->_thermal_solver->_mesh->numberOfNodes() -1]) * initial_power_density;
    MicroCellBoundaryCondition* fixed_flux_boundary_condition = MicroCellBoundaryCondition::getFixedHeatFluxBoundaryConditionFactory(outer_boundary_heat_flux);
    this->_thermal_solver->setBoundaryCondition(fixed_flux_boundary_condition);
    
    
    //Define the kinetics parameters
    this->_kinetics_model = new ReactorKinetics(this,initial_power_density, ReactorKinetics::DelayedPrecursorInitialState::EquilibriumPrecursors);    
    
    //Time stepping parameters
    _power_and_delayed_neutron_record_time_step =  _input_file_reader->getInputFileParameter("Power Record", static_cast<Real>(0.0005) );  //How often to calculate keff and the prompt neutron lifetime
    _kinetics_thermal_sync_time_step = _input_file_reader->getInputFileParameter("Kinetics Thermal Data Sync", static_cast<Real>(20e-6) );      //How often to couple the kinetics and heat transfer routines    
    _end_time = _input_file_reader->getInputFileParameter("Calculation End Time", static_cast<Real>(1.00) );                                    //How many seconds should the simulation last 
        
   

}

void InfiniteCompositeReactor::solveForSteadyStatePowerDistribution(const std::vector<Real> &homogenous_power_density, const Real &initial_power_density)
{
    std::vector<Real> last_power_density(homogenous_power_density);
    Real max_relative_residual, average_residual;
    int initial_power_iteration = 1;

    //Iterate the power distribution and thermal solutions until it 
    do
    {
        //Run a k-effective calculation to get the tally values to get the power density
        this->_monte_carlo_model->getRawKeff();
        std::vector<std::vector<Real>> initial_tally_power_density = this->_monte_carlo_model->getZoneCellRelativePowerDensity();
        std::vector<Real> tally_power_density = this->_thermal_solver->getTallyBasedRepresentativeKernelPowerDistribution(initial_tally_power_density, initial_power_density);
        std::vector<MicroSolution> tally_plot =  this->_thermal_solver->iterateInitialConditions(tally_power_density);
        vector_residuals(last_power_density, tally_power_density,max_relative_residual,average_residual);
        std::cout<< "Max Power Residual: " << max_relative_residual << " Average Power Residual: " << average_residual << std::endl;

        //Save solutions to the graph log
        MicroSolution::plotSolutions(tally_plot,0, this->_results_directory + "initial-solve-" + std::to_string(initial_power_iteration) + ".png"); 
        ++initial_power_iteration;
        last_power_density = tally_power_density;

    } while( (max_relative_residual > 0.005 || average_residual > 0.003) && initial_power_iteration <= 4);
}

void InfiniteCompositeReactor::createOutputFile()
{
    std::ofstream output_file;
    output_file.open( this->_results_directory + this->_data_file, std::ios::out);
    
    output_file << "Iteration,Time [s],Timestep [s],Power [W/m^3],k_eff,k_eff sigma,neutron lifetime [s],Neutron Lifetime sigma [s],Beta_eff,Beta_eff sigma,Run Time [s],Edge Temp [K],Gamma,Power Peaking";
    
    for(size_t index = 1; index <= 6; index++ )
    {
        output_file << ",Group " << index;
    }
    
    output_file << std::endl;    
    output_file.close();
}

/**
 * 
 * @param time
 * @param power
 * @param k_eff
 * @param k_eff_sigma
 * @param neutron_lifetime
 * @param neutron_lifetime_sigma
 * @param beta_eff
 * @param beta_eff_sigma
 * @param hot_temperature
 */
void InfiniteCompositeReactor::saveCurrentData(const Real &time, const Real &power, const Real &k_eff, const Real &k_eff_sigma, const Real &neutron_lifetime, const Real &neutron_lifetime_sigma, const Real &beta_eff, const Real &beta_eff_sigma, const Real &hot_temperature, const Real &gamma)
{
    std::ofstream output_file;
    std::string file_path = this->_results_directory + this->_data_file;
    
    if(!file_exists(file_path))
    {
        this->createOutputFile();
    }
    
    output_file.open( file_path, std::ios::app);
    
    //time elapsed since the simulation started
    std::time_t elapsed_time_since_start =  std::time(nullptr) - InfiniteCompositeReactor::_simulation_start_time;   
    
    Real power_peaking = -1;
    
    if(_monte_carlo_model->_tally_cells)
    {
        std::vector<std::vector<Real>> zones = this->_monte_carlo_model->getZoneCellRelativePowerDensity();
        Real max = zones[0][zones[0].size() - 1];
        Real min = zones[0][0];
        power_peaking = min/max;
    }
    
    output_file << _monte_carlo_number_iterations << "," << time << "," << _monte_carlo_time_iteration << "," << power << "," <<  k_eff << "," << k_eff_sigma << "," << neutron_lifetime << "," << neutron_lifetime_sigma << "," << beta_eff << "," << beta_eff_sigma << "," << elapsed_time_since_start << "," << hot_temperature << "," << gamma << "," << power_peaking;
    
    auto delayed_precursors = _kinetics_model->_delayed_precursors;
    
    for(size_t index = 0; index < delayed_precursors.size(); index++ )
    {
        output_file << "," << delayed_precursors[index];
    }

    output_file << std::endl;    
    output_file.close();
}

InfiniteCompositeReactor::MoneCarloRecalculation InfiniteCompositeReactor::getRecalculationType(const std::string &string_type)
{
    if(string_type == "Temperature")
    {
        return InfiniteCompositeReactor::MoneCarloRecalculation::Temperature;
    }
    else if( string_type == "Time" )
    {
        return InfiniteCompositeReactor::MoneCarloRecalculation::Time;
    }
}