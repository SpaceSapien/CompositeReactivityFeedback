/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousWorth.cpp
 * Author: chris
 * 
 * Created on March 15, 2017, 7:09 AM
 */
#include "InputDataFunctions.h"
#include "HomogenousWorth.h"
#include "InfiniteHomogenousReactor.h"
#include "HomogenousMicroCell.h"
#include "HomogenousMonteCarlo.h"

HomogenousWorth::HomogenousWorth(InfiniteHomogenousReactor* reactor) : WorthStudy(reactor) 
{
    _reactor = reactor;
}

HomogenousWorth::~HomogenousWorth() {}


void HomogenousWorth::setThermalSolver(HomogenousMicroCell* solver)
{
    this->_thermal_solver = solver;
    this->WorthStudy::setThermalSolver(solver);
}

void HomogenousWorth::setMonteCarloModel(HomogenousMonteCarlo* monte_carlo_model)
{
    _monte_carlo_model = monte_carlo_model;
    this->WorthStudy::setMonteCarloModel(monte_carlo_model);
}


void HomogenousWorth::startStudy(const std::string& output_file_name)
{
    //Default geometry sizes
    Real default_sphere_outer_radius = 2e-3;  //meters
    Real default_fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { default_fuel_kernel_outer_radius, default_sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    //Grab geometry from input file
    std::vector<Dimension> dimensions = _reactor->_input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =  _reactor->_input_file_reader->getInputFileParameter("Materials", default_materials);
    _micro_sphere_geometry = new MaterialLibrary::MicroGeometry(materials, dimensions);  
    
    
    std::pair<Real,Real> default_worth_temperature_range = { 300, 4000 };
    std::pair<Real,Real> worth_temperature_range= _reactor->_input_file_reader->getInputFileParameter("Worth Range", default_worth_temperature_range);
    Real starting_temperature = worth_temperature_range.first;
    Real ending_temperature = worth_temperature_range.second;    
    int default_divisions = 12;
    int divisions =  _reactor->_input_file_reader->getInputFileParameter("Worth Temperature Divisions", default_divisions);
    
    Real temperature_increment = (ending_temperature - starting_temperature)/(divisions-1);
    Real current_temperature = starting_temperature;   
    
    _reactor->setMicroSphereGeometry(_micro_sphere_geometry);
    
    this->setMonteCarloModel(new HomogenousMonteCarlo(_reactor, _reactor->_results_directory + "run/")); 
    _monte_carlo_model->_tally_cells = true;
    
    for( int division_counter = 0; division_counter < divisions; ++division_counter)
    {
        std::cout << current_temperature << std::endl;
        this->setThermalSolver(new HomogenousMicroCell(this->_reactor, current_temperature ));
        _reactor->setThermalSolver(_thermal_solver);       
        SimulationResults results = _monte_carlo_model->getRawKeff();          
        this->logData(results, output_file_name);
        
        delete _thermal_solver;        
        _reactor->setThermalSolver(nullptr);
        
        
        current_temperature += temperature_increment;
    }
    
    delete _micro_sphere_geometry;
    _reactor->setMicroSphereGeometry(nullptr);
    delete _monte_carlo_model;
}


