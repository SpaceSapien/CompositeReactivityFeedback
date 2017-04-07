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

#include "InfiniteCompositeReactor.h"
#include "ReactivityInsertion.h"
#include "CompositeMicroCell.h"
#include "InfiniteCompositeWorth.h"
#include "ReactorKinetics.h"
/*
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InputFileParser.h"
#include "CompositeMicroCell.h"*/

InfiniteCompositeReactor::InfiniteCompositeReactor(const std::string &input_file_name ) : Reactor(input_file_name){}

//InfiniteCompositeReactor::InfiniteCompositeReactor(const std::string old_results_folder, Real new_end_time) : Reactor(old_results_folder,new_end_time) {}


void InfiniteCompositeReactor::worthStudy()
{
    std::cout<<"Starting Worth Study"<<std::endl;
    InfiniteCompositeWorth* study = new InfiniteCompositeWorth(this);
    
    //Setting certain times for the transient time to track the tallies
    this->_transient_time = -4.0;
    study->startStudy(true,false,"fuel-worth.csv");
    this->_transient_time = -3.0;
    study->startStudy(false,true,"moderator-worth.csv");
    this->_transient_time = -2.0;
    study->startStudy(true,true,"worth.csv");
    
    delete study;
    
}

/**
 * Deconstructor - self explanitory
 */
InfiniteCompositeReactor::~InfiniteCompositeReactor(){}

/**
 * Create objects, solve for initial conditions, and otherwise do all work
 * related to solving for pre-transient conditions
 */
void InfiniteCompositeReactor::initializeReactorProblem()
{
    //Call the parent constructor
    Reactor::initializeReactorProblem(); 
          
    //Define the heat transfer settings
    Real initial_power_density =  _input_file_reader->getInputFileParameter("Starting Power Density",static_cast<Real>(200e6) ); // W/m^3 averaged over the entire micro sphere
    Real initial_outer_shell_temperature = _input_file_reader->getInputFileParameter("Kernel Outer Temperature",800); // Kelvin
    
    //Set up a boundary condition 
    MicroCellBoundaryCondition* fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(initial_outer_shell_temperature);
    //Then crate the MicoCell thermal solver, set the boundary condition and then iterate steady state condition 
    this->setThermalSolver(new CompositeMicroCell(this, initial_outer_shell_temperature));
    this->_thermal_solver->setBoundaryCondition(fixed_temperature_boundary_condition);
    
    std::vector<Real> homogenous_power_density = this->_thermal_solver->getRespresentativePowerDistribution(initial_power_density);
    std::vector<MicroSolution> plot =  this->_thermal_solver->iterateInitialConditions(homogenous_power_density);
    MicroSolution::plotSolutions(plot,0, this->_results_directory + "initial-pre-tally-solve.png");
    
    //Define the Monte Carlo Parameters
    this->setMoteCarloModel(new InfiniteCompositeMonteCarlo(this, this->_results_directory + "run/"));   
    
    if(_monte_carlo_model->_tally_cells)
    {
        this->solveForSteadyStatePowerDistribution(homogenous_power_density,initial_power_density);            
    }    
        
    this->_monte_carlo_model->updateAdjustedCriticalityParameters();
    
    
    //Define the kinetics parameters
    this->setKineticsModel(new ReactorKinetics(this,initial_power_density, ReactorKinetics::DelayedPrecursorInitialState::EquilibriumPrecursors, this->_monte_carlo_model->_current_beta_eff));    
    
    //Reset the time to zero and then grab the current solution for the Problem
    _thermal_solver->_current_time = 0;    
     //Add the starting MicroSolution
    _plot_solutions.push_back( _thermal_solver->getCurrentMicrosolution() );
        
    
    //Setting Up the Transient BC    
    std::string transient_boundary_condition = _input_file_reader->getInputFileParameter(std::string("Transient BC"), std::string("FixedHeatFlux") );
    MicroCellBoundaryCondition::BoundaryType transient_boundary_type = MicroCellBoundaryCondition::getBoudaryTypeFromString(transient_boundary_condition);
    
    if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::FixedHeatFlux )
    {
        //Get the steady state heat flux and set it as the boundary condition
        Real outer_boundary_heat_flux = sphere_volume(this->_thermal_solver->_mesh->_outer_radius[this->_thermal_solver->_mesh->numberOfNodes() -1]) * initial_power_density;
        MicroCellBoundaryCondition* fixed_flux_boundary_condition = MicroCellBoundaryCondition::getFixedHeatFluxBoundaryConditionFactory(outer_boundary_heat_flux);
        this->_thermal_solver->setBoundaryCondition(fixed_flux_boundary_condition);
    }
    else if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::FixedTemperature )
    {
        //do nothing as this BC is already set
    }
    else if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::ReflectedHeatFlux )
    {
        MicroCellBoundaryCondition* reflected_boundary = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();
        this->_thermal_solver->setBoundaryCondition(reflected_boundary);
    }
    else
    {
        throw "Unknown BC for transient";
    }
    
    //Define the transient parameters
    this->setReactivityInsertionModel( new ReactivityInsertion(this) );  
}

/**
 * Solve for the power density and the steady state temperature for that power density. Then feed that into a MC run to get
 * another power density. If the power density has converged then we have the proper power density
 * @param homogenous_power_density
 * @param initial_power_density
 */
void InfiniteCompositeReactor::solveForSteadyStatePowerDistribution(const std::vector<Real> &homogenous_power_density, const Real &initial_power_density)
{
    std::vector<Real> last_power_density(homogenous_power_density);
    Real max_relative_residual, average_residual;
    int initial_power_iteration = 1;
    
     Real power_residual = _input_file_reader->getInputFileParameter("Initial Solve Power Residual", static_cast<Real>(0.005) );

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

    } while( (max_relative_residual > power_residual || average_residual > (power_residual /2) ) && initial_power_iteration <= 4);
}

void InfiniteCompositeReactor::setThermalSolver(CompositeMicroCell* solver)
{
    this->_thermal_solver = solver;
    this->Reactor::setThermalSolver(solver);
}

 void InfiniteCompositeReactor::setMoteCarloModel(InfiniteCompositeMonteCarlo* mc_model)
 {
     this->_monte_carlo_model = mc_model;
     this->Reactor::setMoteCarloModel(mc_model);
 }