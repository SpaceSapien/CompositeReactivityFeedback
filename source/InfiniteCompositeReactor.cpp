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
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include <math.h>
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


InfiniteCompositeReactor::InfiniteCompositeReactor() 
{    
    //Create the Results Folder
    time_t run_identification_number = std::time(nullptr);
    _results_directory =  "results/" + std::to_string(run_identification_number) + "/";    
    std::string folder_command = "mkdir -p " + _results_directory;
    system( folder_command.c_str());
    
    initializeInifiniteCompositeReactorProblem();
}


void InfiniteCompositeReactor::simulate()
{
    //Save the initial conditions in the plot_solutions array
    MicroSolution solution =  this->_thermal_solver->getInitialConditions();
    _plot_solutions.push_back(solution);
    
    Real inner_time_step = 0;
    
    //Simulate the transient the outer loop is the monte carlo simulation
    for( Real transient_time = 0; transient_time < _end_time; transient_time+= inner_time_step)
    {
        //Gather the parameters from the monte carlo model 
        //The Monte Carlo model is run on the outer loop
        Real prompt_removal_lifetime = _monte_carlo_model->_current_prompt_neutron_lifetime;
        Real k_eff = _monte_carlo_model->_current_k_eff; 
        Real lambda = prompt_removal_lifetime/k_eff;
        std::vector<std::pair<FissionableIsotope,Real> > fission_listing = _monte_carlo_model->_fission_tally_listing;
        Real reactivity = 10000.0*(k_eff - 1.0)/k_eff;
        
        std::pair<Real,Real> reactivity_pair = { transient_time , reactivity };
        _reactivity_record.push_back( reactivity_pair );
        std::pair<Real,Real> prompt_removal_lifetime_pair = { transient_time, prompt_removal_lifetime };
        _prompt_life_time_record.push_back(prompt_removal_lifetime_pair);
       
         
        for( inner_time_step = 0 ; inner_time_step < _monte_carlo_time_iteration ; inner_time_step += _kinetics_time_iteration)
        {
            //Solve the kinetics model
            Real current_power = _kinetics_model->solveForPower(_kinetics_time_iteration, k_eff,lambda,fission_listing, _power_record, _delayed_record);
            std::vector<Real> power_distribition = getPowerDistribution(_thermal_solver->_solver_settings._radial_mesh,_micro_sphere_geometry->getFuelKernelRadius(), current_power);

            //Get the thermal solution
            solution =  _thermal_solver->solve( _kinetics_time_iteration, power_distribition); 
            //Add the current time steps solution to the mix
        }
        
        solution.plot( this->_results_directory + "solution-" + std::to_string(_thermal_solver->_current_time)  + ".png");
        _plot_solutions.push_back(solution); 
        
        //if there is still enough time left to do another monte carlo time iteration
        if(transient_time + inner_time_step < _end_time)
        {
            _monte_carlo_model->updateAdjustedCriticalityParameters();
        }
         
        
    }
    
    MicroSolution::saveSolutions( _plot_solutions, this->_results_directory );
    MicroSolution::plotSolutions( _plot_solutions,8 , this->_results_directory + "solutions-graph.png");
    PythonPlot::plotData( _power_record, "Time [s]", "Power Density [W/m^3]","","Power vs. Time", this->_results_directory + "power-graph.png");
    PythonPlot::plotData( _prompt_life_time_record, "Time [s]", "Prompt Neutron Lifetime [s]","","Prompt Neutron Lifetime vs. Time", this->_results_directory + "prompt-neutron-lifetime-graph.png");
    PythonPlot::plotData( _reactivity_record, "Time [s]", "Excess Reactivity [pcm]","","Excess Reactivity vs. Time", this->_results_directory + "excess-reactivity-graph.png");
    PythonPlot::plotData( _delayed_record, "Time [s]", "Delayed Precursors", {} , "Keff vs. Delayed Precursors", "delayed-precursors.png");
}

InfiniteCompositeReactor::~InfiniteCompositeReactor()
{
    delete _micro_sphere_geometry;
    delete _thermal_solver;
    delete _kinetics_model;
    delete _monte_carlo_model;
}
    
void InfiniteCompositeReactor::initializeInifiniteCompositeReactorProblem()
{
    //Define our geometry
    Real sphere_outer_radius = 2e-3;  //meters
    Real fuel_kernel_outer_radius = 4e-4;
    std::vector<Dimension> dimensions = { fuel_kernel_outer_radius, sphere_outer_radius };
    std::vector<Materials> materials = { Materials::UO2, Materials::C }; 
    this->_micro_sphere_geometry = new MicroGeometry(materials, dimensions);    
    
    //Define the heat transfer settings
    Real initial_power_density = 200e6; // W/m^3 averaged over the entire micro sphere
    Real initial_outer_shell_temperature = 800;
    MicroCellBoundaryCondition fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedBoundaryCondition(initial_outer_shell_temperature);
    this->_thermal_solver = new MicroCell(*this->_micro_sphere_geometry,initial_outer_shell_temperature);
    this->_thermal_solver->setBoundaryCondition(fixed_temperature_boundary_condition);
    std::vector<MicroSolution> plot =  this->_thermal_solver->iterateInitialConditions(initial_power_density);
    MicroSolution::plotSolutions(plot,0, this->_results_directory + "initial-solve.png");    
    Real outer_boundary_spatial_derivate = this->_thermal_solver->getOuterDerivative();
    MicroCellBoundaryCondition reflected_boundary_condition = MicroCellBoundaryCondition::getFixedDerivativeBoundaryCondition(outer_boundary_spatial_derivate);
    //MicroCellBoundaryCondition reflected_boundary_condition = MicroCellBoundaryCondition::getRefelectedBoundaryCondition();
    this->_thermal_solver->setBoundaryCondition(reflected_boundary_condition);
    
    //Define the Monte Carlo Parameters
    Real starting_k_eff = 1.01;    
    this->_monte_carlo_model = new ReactorMonteCarlo(this->_thermal_solver, starting_k_eff, this->_results_directory + "run/");   
    
    //Define the kinetics parameters
    std::vector<std::pair<FissionableIsotope,Real> > fission_listing = this->_monte_carlo_model->_fission_tally_listing;
    this->_kinetics_model = new ReactorKinetics(initial_power_density, ReactorKinetics::DelayedPrecursorInitialState::EquilibriumPrecursors);
    
    
    //Time stepping parameters
    _monte_carlo_time_iteration = 0.01;  //How often to calculate keff and the prompt neutron lifetime
    _kinetics_time_iteration = 0.0002;   //How often to couple the kinetics and heat transfer routines    
    _end_time = 0.05;                    //How many seconds should the simulation last 
    
}