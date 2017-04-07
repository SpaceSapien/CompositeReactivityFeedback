/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinReactor.cpp
 * Author: chris
 * 
 * Created on March 22, 2017, 10:26 PM
 */

#include "FuelPinReactor.h"
#include "ReactivityInsertion.h"
#include "ReactorKinetics.h"
#include <cmath>

FuelPinReactor::FuelPinReactor(const std::string &input_file_name) : Reactor(input_file_name) 
{
    
    
}




FuelPinReactor::DimensionalTreatment FuelPinReactor::getDimensionalTreatment(const std::string input_string)
{
    if(input_string == "HomogenousNeutronicsAndHeatTransfer" )
    {
        return FuelPinReactor::HomogenousNeutronicsAndHeatTransfer;
    }
    else if(input_string == "HomogenousNeutronics" )
    {
        return FuelPinReactor::HomogenousNeutronics;
    }
    else if(input_string == "FullHeterogeneous")
    {
        return FuelPinReactor::FullHeterogeneous;
    }
}



void FuelPinReactor::initializeReactorProblem()
{
    //Call the parent constructor
    Reactor::initializeReactorProblem(); 
    
    //Create a series of macro_cells
    _number_macro_cells = _input_file_reader->getInputFileParameter("Number Macro Cells", 0 );
          
    //Define the heat transfer settings
    Real initial_power_density =  _input_file_reader->getInputFileParameter("Starting Power Density",static_cast<Real>(200e6) ); // W/m^3 averaged over the entire micro sphere
    Real initial_outer_shell_temperature = _input_file_reader->getInputFileParameter("Kernel Outer Temperature",800); // Kelvin
    
    //Set up a boundary condition 
    MicroCellBoundaryCondition* fixed_temperature_boundary_condition = MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(initial_outer_shell_temperature);
    //Then crate the MicoCell thermal solver, set the boundary condition and then iterate steady state condition 
    this->setThermalSolver(new CylindricalMicroCell(this, initial_outer_shell_temperature));
    this->_thermal_solver->setInnerBoundaryCondition(fixed_temperature_boundary_condition);
    
    std::vector<Real> homogenous_power_density = this->_thermal_solver->getRespresentativePowerDistribution(initial_power_density);
    std::vector<MicroSolution> plot =  this->_thermal_solver->iterateInitialConditions(homogenous_power_density);
    MicroSolution::plotSolutions(plot,0, this->_results_directory + "initial-pre-tally-solve.png");
    
    //Define the Monte Carlo Parameters
    this->setMoteCarloModel(new FuelPinMonteCarlo(this, this->_results_directory + "run/"));   
    
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
        Real outer_boundary_heat_flux = M_PI * ( _thermal_solver->_mesh->_outer_radius.back() * _thermal_solver->_mesh->_outer_radius.back() - _thermal_solver->_mesh->_inner_radius.front() * _thermal_solver->_mesh->_inner_radius.front() ) * initial_power_density;
        MicroCellBoundaryCondition* fixed_flux_boundary_condition = MicroCellBoundaryCondition::getFixedHeatFluxBoundaryConditionFactory(outer_boundary_heat_flux);
        this->_thermal_solver->setInnerBoundaryCondition(fixed_flux_boundary_condition);
    }
    else if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::FixedTemperature )
    {
        //do nothing as this BC is already set
    }
    else if(transient_boundary_type == MicroCellBoundaryCondition::BoundaryType::ReflectedHeatFlux )
    {
        MicroCellBoundaryCondition* reflected_boundary = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();
        this->_thermal_solver->setInnerBoundaryCondition(reflected_boundary);
    }
    else
    {
        throw "Unknown BC for transient";
    }
    
    //Define the transient parameters
    this->setReactivityInsertionModel( new ReactivityInsertion(this) );  
    
    std::string dimensionality_string = _input_file_reader->getInputFileParameter("Fuel Pin Dimensionality", std::string("HomogenousNeutronics") );
    _dimensionality = FuelPinReactor::getDimensionalTreatment(dimensionality_string);

}

void FuelPinReactor::worthStudy()
{
    // Not sure if a worth study is relevant?
}

FuelPinReactor::~FuelPinReactor() {}
 

void FuelPinReactor::setThermalSolver(CylindricalMicroCell* solver)
{
    this->_thermal_solver = solver;
    this->Reactor::setThermalSolver(solver);
}

 void FuelPinReactor::setMoteCarloModel(FuelPinMonteCarlo* mc_model)
 {
     this->_monte_carlo_model = mc_model;
     this->Reactor::setMoteCarloModel(mc_model);
 }
 
 void FuelPinReactor::solveForSteadyStatePowerDistribution(const std::vector<Real> &homogenous_power_density, const Real &initial_power_density)
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