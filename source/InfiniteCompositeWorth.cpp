/* 
 * File:   infiniteCompositeWorth.cpp
 * Author: chris
 * 
 * Created on March 15, 2017, 4:42 AM
 */

#include "InputDataFunctions.h"
#include "InfiniteCompositeWorth.h"
#include "InfiniteCompositeReactor.h"
#include "CompositeMicroCell.h"

InfiniteCompositeWorth::InfiniteCompositeWorth(InfiniteCompositeReactor* reactor) : WorthStudy(reactor)
{
    this->_reactor = reactor;
}

InfiniteCompositeWorth::~InfiniteCompositeWorth() {}

void InfiniteCompositeWorth::startStudy(const std::string &output_file_name)
{
    this->startStudy(true,true,output_file_name);
}

void InfiniteCompositeWorth::setThermalSolver(CompositeMicroCell* thermal_solver)
{
    this->_thermal_solver = thermal_solver;
    this->WorthStudy::setThermalSolver(thermal_solver);
}

void InfiniteCompositeWorth::startStudy(const bool &vary_fuel_temperature,const bool &vary_matrix_temperature,const std::string &output_file_name)
{
    //Default geometry sizes
    Real default_sphere_outer_radius = 2e-3;  //meters
    Real default_fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { default_fuel_kernel_outer_radius, default_sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    //Grab geometry from input file
    std::vector<Dimension> dimensions = _reactor->_input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =  _reactor->_input_file_reader->getInputFileParameter("Materials", default_materials);
    
    std::pair<Real,Real> default_worth_temperature_range = { 300, 4000 };
    std::pair<Real,Real> worth_temperature_range= _reactor->_input_file_reader->getInputFileParameter("Worth Range", default_worth_temperature_range);
    Real starting_temperature = worth_temperature_range.first;
    Real ending_temperature = worth_temperature_range.second;
    
    int default_divisions = 12;
    int divisions =  _reactor->_input_file_reader->getInputFileParameter("Worth Temperature Divisions", default_divisions);
    
    Real temperature_increment = (ending_temperature - starting_temperature)/(divisions-1);
    Real current_temperature = starting_temperature;
    
    _micro_sphere_geometry = new MaterialLibrary::MicroGeometry(materials, dimensions);  
    _reactor->setMicroSphereGeometry(_micro_sphere_geometry);
    _monte_carlo_model = new InfiniteCompositeMonteCarlo(_reactor, _reactor->_results_directory + "run/"); 
    _monte_carlo_model->_tally_cells = true;
    
    for( int division_counter = 0; division_counter < divisions; ++division_counter)
    {
        std::cout << current_temperature << std::endl;
        
        if( vary_fuel_temperature)
        {
            this->setThermalSolver(new CompositeMicroCell(_reactor, current_temperature ));
        }
        else
        {
            this->setThermalSolver(new CompositeMicroCell(_reactor,starting_temperature ));
        }
        
        if( vary_matrix_temperature )
        {
            _thermal_solver->setOuterMaterialTemperature(current_temperature);
        }
        else
        {
            _thermal_solver->setOuterMaterialTemperature(starting_temperature);
        }
        
        
        
        _reactor->setThermalSolver(_thermal_solver);
        
        
        SimulationResults results = _monte_carlo_model->getRawKeff();          
        this->logData(results, output_file_name);
        
        
        delete _thermal_solver;
        
        _reactor->_thermal_solver = nullptr;
        _reactor->Reactor::_thermal_solver = nullptr;
        
        current_temperature += temperature_increment;
    }
    
    delete _micro_sphere_geometry;
    _reactor->_micro_sphere_geometry = nullptr;
    delete _monte_carlo_model;
}





