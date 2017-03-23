/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InfiniteHomogenousReactor.cpp
 * Author: chris
 * 
 * Created on March 13, 2017, 11:54 AM
 */

#include "InfiniteHomogenousReactor.h"
#include "HomogenousWorth.h"
#include "ReactivityInsertion.h"

InfiniteHomogenousReactor::InfiniteHomogenousReactor(const std::string &input_file_name) : Reactor(input_file_name) {}

InfiniteHomogenousReactor::~InfiniteHomogenousReactor() {}




void InfiniteHomogenousReactor::initializeReactorProblem()
{
    //Call the parent constructor
    Reactor::initializeReactorProblem(); 
          
    //Define the heat transfer settings
    Real initial_power_density =  _input_file_reader->getInputFileParameter("Starting Power Density",static_cast<Real>(200e6) ); // W/m^3 averaged over the entire micro sphere
    Real initial_temperature = _input_file_reader->getInputFileParameter("Kernel Outer Temperature",800); // Kelvin
    
    this->setThermalSolver(new HomogenousMicroCell(this, initial_temperature));
    
    //Define the Monte Carlo Parameters
    this->setMoteCarloModel(new HomogenousMonteCarlo(this, this->_results_directory + "run/"));   
    this->_monte_carlo_model->updateAdjustedCriticalityParameters();
    
    //Setting Up the Transient BC    
    std::string transient_boundary_condition = _input_file_reader->getInputFileParameter(std::string("Transient BC"), std::string("FixedHeatFlux") );
    MicroCellBoundaryCondition::BoundaryType transient_boundary_type = MicroCellBoundaryCondition::getBoudaryTypeFromString(transient_boundary_condition);
    
    if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::FixedHeatFlux )
    {
        //Get the steady state heat flux and set it as the boundary condition
        Real outer_boundary_heat_flux = initial_power_density;
        MicroCellBoundaryCondition* fixed_flux_boundary_condition = MicroCellBoundaryCondition::getFixedHeatFluxBoundaryConditionFactory(outer_boundary_heat_flux);
        this->_thermal_solver->setBoundaryCondition(fixed_flux_boundary_condition);
    }
    else if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::FixedTemperature )
    {
        MicroCellBoundaryCondition* fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(initial_temperature);
        this->_thermal_solver->setBoundaryCondition(fixed_temperature_boundary_condition);
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
    
    this->_thermal_solver->_outward_current_power = initial_power_density;
}
void InfiniteHomogenousReactor::worthStudy()
{
    HomogenousWorth* worth_study = new HomogenousWorth(this);
    
    worth_study->startStudy();
    
    delete worth_study;
}
/*void InfiniteHomogenousReactor::postSimulationProcessing()
{
    
}
void InfiniteHomogenousReactor::monteCarloTimeStepSimulationProcessing()
{
    
}*/

void InfiniteHomogenousReactor::setThermalSolver(HomogenousMicroCell* solver)
{
    this->_thermal_solver = solver;
    this->Reactor::setThermalSolver(solver);
}

