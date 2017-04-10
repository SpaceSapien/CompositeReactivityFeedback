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
#include "EnumsAndFunctions.h"
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
    
    // Set the dimensionality
    std::string dimensionality_string = _input_file_reader->getInputFileParameter("Fuel Pin Dimensionality", std::string("HomogenousNeutronics") );
    _dimensionality = FuelPinReactor::getDimensionalTreatment(dimensionality_string);
    
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
    
    

}



bool FuelPinReactor::significantTemperatureDifference(MicroSolution* comparison)
{
    bool significan_difference_in_macrocell = this->Reactor::significantTemperatureDifference(comparison);
    
    if(significan_difference_in_macrocell)
    {
        return true;
    }    
        
    
    //Add in the cell data
    switch( _dimensionality )
    {
        //Heterogeneous MicroCells
        case FuelPinReactor::HomogenousNeutronics :
        case FuelPinReactor::FullHeterogeneous :
        {
            std::vector<Real> run_data = getMicroCellTemperatures();
            
            if(this->significantTemperatureDifferenceInMicroCell(run_data))
            {
                return true;
            }
            
            break;
        }
    }
    
    return false;
}

/**
 * Comparing the current microsolution to the one on record as the last k-eignvalue time step microsolution
 * 
 * @param comparison the current microsolution
 * @return true or false if there is a significant enought temperature change to warrant a new k-eignevalue calc
 */
bool FuelPinReactor::significantTemperatureDifferenceInMicroCell(std::vector<Real> comparison)
{
     std::vector<Real> last_solution_microcell_temperatures = _microsolver_average_temperatures.back();
     
     Real max_difference;
     Real average_difference;
     
     vector_difference(comparison, last_solution_microcell_temperatures, max_difference, average_difference);
     
     if( max_difference > _maximum_allowed_temperature_difference || average_difference > _maximum_allowed_average_temperature_difference )
     {
         return true;
     }
     else
     {
         return false;
     }     
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
        
        for(std::size_t index = 0; index < _thermal_solver->_micro_scale_solvers.size(); ++index)
        {
            delete _thermal_solver->_micro_scale_solvers[index];
            _thermal_solver->_micro_scale_solvers[index] = nullptr;            
        }
        
        std::vector<MicroSolution> tally_plot =  this->_thermal_solver->iterateInitialConditions(tally_power_density);
        vector_residuals(last_power_density, tally_power_density,max_relative_residual,average_residual);
        std::cout<< "Max Power Residual: " << max_relative_residual << " Average Power Residual: " << average_residual << std::endl;

        //Save solutions to the graph log
        MicroSolution::plotSolutions(tally_plot,0, this->_results_directory + "initial-solve-" + std::to_string(initial_power_iteration) + ".png"); 
        ++initial_power_iteration;
        last_power_density = tally_power_density;

    } while( (max_relative_residual > power_residual || average_residual > (power_residual /2) ) && initial_power_iteration <= 4);
    
    
    
}
 
 
void FuelPinReactor::monteCarloTimeStepSimulationDataProcessing()
{
    this->Reactor::monteCarloTimeStepSimulationDataProcessing();
     
    //Add in the macro cell temperature change data
    switch( _dimensionality )
    {
        //Heterogeneous MicroCells
        case FuelPinReactor::HomogenousNeutronics :
        case FuelPinReactor::FullHeterogeneous :
        {
            std::vector<Real> run_data = getMicroCellTemperatures();
            _microsolver_average_temperatures.push_back(run_data);
            this->saveMicroScaleDataToFile();
            
            break;
        }
    }
}
 
//Keep a record of the micro cell temperatures
 std::vector<Real> FuelPinReactor::getMicroCellTemperatures()
 {
     std::vector<Real> run_data;
     
    for( auto cell : this->_thermal_solver->_micro_scale_solvers )
    {
        Real temperature, volume;
        cell->getAverageZoneTemperature(1,temperature,volume);
        run_data.push_back(temperature);
    }
     
     return run_data;
 }
 
void FuelPinReactor::saveMicroScaleDataToFile()
{
    std::ofstream output_file;
    std::string output_file_path = _results_directory + "microscale-aggregate-data.csv";
        
    if(!file_exists(output_file_path))
    {
        
        
        output_file.open( output_file_path, std::ios::app);

        output_file << "Time [s]";
        
        for( int index = 0; index < _thermal_solver->_micro_scale_solvers.size(); ++index )
        {
            output_file << ",Position-" << index << " [m],Temperature-"<< index << " [K],Integrated-Power-" << index << " [W-s],Current-Power-" << index <<" [W]";
        }
    
        output_file << std::endl;    
    }
    else
    {
        output_file.open( output_file_path, std::ios::app);
    }
    
    output_file << _transient_time;
    
        
    for( auto cell : this->_thermal_solver->_micro_scale_solvers )
    {
        Real temperature,volume;
        
        cell->getAverageZoneTemperature(1,temperature,volume);
        
        output_file << "," << cell->_macro_scale_position;
        output_file << "," << temperature;
        output_file << "," << cell->_outward_integrated_power;
        output_file << "," << cell->_outward_current_power;       
    }
    output_file << std::endl;  
    output_file.close();
    
    
    
    //Fine data
    output_file_path = _results_directory + "microscale-fine-data.csv";
        
    if(!file_exists(output_file_path))
    {
        
        
        output_file.open( output_file_path, std::ios::app);
        output_file << "Time [s],Cell,";
        
        for( int index = 0; index < _thermal_solver->_micro_scale_solvers[0]->_mesh->_position.size(); ++index )
        {
            output_file << ",Radius-" << index << " [m],Temperature-"<< index << " [K]";
        }
    
        output_file << std::endl;    
    }
    else
    {
        output_file.open( output_file_path, std::ios::app);
    }
    
    
    
    
    for( auto cell : _thermal_solver->_micro_scale_solvers )
    {
        output_file << _transient_time;
        
        for(int index = 0; index < cell->_mesh->_position.size(); ++index)
        {
            output_file << "," << cell->_mesh->_position[index] << "," << cell->_solution[index];
        }
        
        output_file << std::endl;  
    }
    
    output_file.close();
 }