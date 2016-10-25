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
#include "ExplicitSolverSettings.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InputFileParser.h"

const std::time_t InfiniteCompositeReactor::_simulation_start_time = std::time(nullptr);;

InfiniteCompositeReactor::InfiniteCompositeReactor(const std::string &input_file_name ) 
{   
    
    
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
    
    
    
    
    
    createOutputFile();
    initializeInifiniteCompositeReactorProblem();
}


void InfiniteCompositeReactor::simulate()
{
   
    //Simulate the transient the outer loop is the monte carlo simulation
    for( this->_transient_time = 0; _transient_time < _end_time; _transient_time += _inner_time_step)
    {
        this->timeIterationInnerLoop();
        this->monteCarloTimeStepSimulationProcessing();
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
    Real reactivity_uncertainty = std::sqrt( 2 ) * k_eff_sigma;
    Real reactivity_pcm = 10000.0 * reactivity;
    Real reactivity_pcm_uncertainty = 10000.0 * reactivity_uncertainty; 
    Real reactivity_cents = reactivity / beta_eff;
    Real reactivity_cents_uncertainty = std::sqrt( (reactivity_uncertainty / reactivity) * (reactivity_uncertainty / reactivity) + ( beta_eff_sigma / beta_eff) * ( beta_eff_sigma / beta_eff) );

    std::tuple<Real,Real,Real> reactivity_pcm_pair = std::make_tuple(_transient_time, reactivity_pcm, reactivity_pcm_uncertainty);
    _reactivity_pcm_record.push_back(reactivity_pcm_pair);        

    std::tuple<Real,Real,Real> reactivity_cents_pair = std::make_tuple(_transient_time, reactivity_cents, reactivity_cents_uncertainty);
   _reactivity_cents_record.push_back(reactivity_cents_pair);

    Real hottest_temperature = this->_thermal_solver->_solution[0]; 

    this->saveCurrentData(_transient_time, current_power, k_eff, k_eff_sigma, lambda, lambda_sigma, beta_eff, beta_eff_sigma, hottest_temperature);
    
    MicroSolution solution = this->_thermal_solver->getCurrentMicrosolution();
    
    solution.plot( this->_results_directory + "solution-" + std::to_string( _transient_time )  + ".png", 800, 3000);
    _plot_solutions.push_back(solution); 
        
}

void InfiniteCompositeReactor::postSimulationProcessing()
{
    //Post Processing graph creation
    MicroSolution::saveSolutions( _plot_solutions, this->_results_directory );
    MicroSolution::plotSolutions( _plot_solutions, 4 , this->_results_directory + "solutions-graph.png");
    PythonPlot::plotData(      _power_record,            "Time [s]", "Power Density [W/m^3]",      "", "Power vs. Time",                   this->_results_directory + "power-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _prompt_life_time_record, "Time [s]", "Prompt Neutron Lifetime [s]","", "Prompt Neutron Lifetime vs. Time", this->_results_directory + "prompt-neutron-lifetime-graph.png", {0, _end_time} );
    PythonErrorPlot::plotData( _k_eff_record,            "Time [s]", "K effective",    "",             "K-eff vs. Time",                   this->_results_directory + "k-eff-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_cents_record, "Time [s]", "Reactivity [Dollars]",    "",    "Reactivity vs. Time",              this->_results_directory + "reactivity-cents-graph.png",        {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_pcm_record,   "Time [s]", "Reactivity [pcm]",    "",        "Reactivity vs. Time",              this->_results_directory + "reactivity-pcm-graph.png",          {0, _end_time} );
    PythonErrorPlot::plotData( _beta_eff_record,         "Time [s]", "Beta effective",    "",          "Beta-eff vs. Time",                this->_results_directory + "beta-eff-graph.png",                {0, _end_time} );
    PythonPlot::plotData(      _delayed_record,          "Time [s]", "Delayed Precursors",         {}, "Keff vs. Delayed Precursors",      this->_results_directory + "delayed-precursors.png",            {0, _end_time} );
    PythonPlot::createPlots();
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
}

/*InfiniteCompositeReactor::temperatureIterationInnerLoop()
{
    
}*/


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
    //This command establishes the logging of python plot commands. Useful for debugging
    PythonPlot::_log_file = this->_results_directory + "graph_log.log";
    
    Real sphere_outer_radius, fuel_kernel_outer_radius; 
    
    //Define our geometry
    sphere_outer_radius = 2e-3;  //meters
    fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { fuel_kernel_outer_radius, sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    std::vector<Dimension> dimensions = _input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =   _input_file_reader->getInputFileParameter("Materials", default_materials);
    
    this->_micro_sphere_geometry = new MicroGeometry(materials, dimensions);    
    
    
    
    //Define the heat transfer settings
    Real initial_power_density =  _input_file_reader->getInputFileParameter("Starting Power Density",static_cast<Real>(200e6) ); // W/m^3 averaged over the entire micro sphere
    Real initial_outer_shell_temperature = 800;//_input_file.getInputFileParameter("Kernel Outer Temperature",800); // Kelvin
    
    MicroCellBoundaryCondition fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedBoundaryCondition(initial_outer_shell_temperature);
    
    this->_thermal_solver = new MicroCell(this, initial_outer_shell_temperature);
    this->_thermal_solver->setBoundaryCondition(fixed_temperature_boundary_condition);
    std::vector<MicroSolution> plot =  this->_thermal_solver->iterateInitialConditions(initial_power_density);
    MicroSolution::plotSolutions(plot,0, this->_results_directory + "initial-solve.png");    
    Real outer_boundary_spatial_derivate = this->_thermal_solver->getOuterDerivative();
    MicroCellBoundaryCondition reflected_boundary_condition = MicroCellBoundaryCondition::getFixedDerivativeBoundaryCondition(outer_boundary_spatial_derivate);
    //MicroCellBoundaryCondition reflected_boundary_condition = MicroCellBoundaryCondition::getRefelectedBoundaryCondition();
    this->_thermal_solver->setBoundaryCondition(reflected_boundary_condition);
    
    //Define the Monte Carlo Parameters
    Real starting_k_eff =  _input_file_reader->getInputFileParameter("Starting K-eff",static_cast<Real>(1.01) );    
    this->_monte_carlo_model = new ReactorMonteCarlo(this, starting_k_eff, this->_results_directory + "run/");   
    
    //Define the kinetics parameters
    this->_kinetics_model = new ReactorKinetics(this,initial_power_density, ReactorKinetics::DelayedPrecursorInitialState::EquilibriumPrecursors);    
    
    //Time stepping parameters
    _power_and_delayed_neutron_record_time_step =  _input_file_reader->getInputFileParameter("Power Record", static_cast<Real>(0.0005) );  //How often to calculate keff and the prompt neutron lifetime
    _monte_carlo_time_iteration =  _input_file_reader->getInputFileParameter("Monte Carlo Recalculation Timestep", static_cast<Real>(0.01) );  //How often to calculate keff and the prompt neutron lifetime
    _kinetics_thermal_sync_time_step = _input_file_reader->getInputFileParameter("Kinetics Thermal Data Sync", static_cast<Real>(20e-6) );      //How often to couple the kinetics and heat transfer routines    
    _end_time = _input_file_reader->getInputFileParameter("Calculation End Time", static_cast<Real>(1.00) );                                    //How many seconds should the simulation last 
    
    //Resetting the timer to zero so that the timer reads t = 0 when the transient starts
    _kinetics_model->_current_time = 0;
}

void InfiniteCompositeReactor::createOutputFile()
{
    std::ofstream output_file;
    output_file.open( this->_results_directory + this->_data_file, std::ios::out);
    
    output_file << "Time [s], Power [W/m^3], k_eff, k_eff sigma, neutron lifetime [s], Neutron Lifetime sigma [s], Beta_eff, Beta_eff sigma, Run Time [s], Edge Temp [K]";
    
    for(size_t index = 1; index <= 6; index++ )
    {
        output_file << ", Group " << index;
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
void InfiniteCompositeReactor::saveCurrentData(const Real &time, const Real &power, const Real &k_eff, const Real &k_eff_sigma, const Real &neutron_lifetime, const Real &neutron_lifetime_sigma, const Real &beta_eff, const Real &beta_eff_sigma, const Real &hot_temperature)
{
    std::ofstream output_file;
    output_file.open( this->_results_directory + this->_data_file, std::ios::app);
    
    //time elapsed since the simulation started
    std::time_t elapsed_time_since_start =  std::time(nullptr) - InfiniteCompositeReactor::_simulation_start_time;
    
    
    
    output_file << time << ", " << power << ", " <<  k_eff << ", " << k_eff_sigma << ", " << neutron_lifetime << ", " << neutron_lifetime_sigma << ", " << beta_eff << ", " << beta_eff_sigma << ", " << elapsed_time_since_start << ", " << hot_temperature;
    
    size_t delayed_size = _delayed_record.size();
    
    if(delayed_size > 0 )
    {
    
        auto delayed_precursors = _delayed_record[delayed_size -1].second;

        for(size_t index = 0; index < delayed_precursors.size(); index++ )
        {
            output_file << ", " << delayed_precursors[index];
        }
    }
    output_file << std::endl;    
    output_file.close();
}