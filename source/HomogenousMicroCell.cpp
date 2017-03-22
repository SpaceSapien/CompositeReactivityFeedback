/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMicroCell.cpp
 * Author: chris
 * 
 * Created on March 15, 2017, 7:05 AM
 */
#include <iostream>
#include "HomogenousMicroCell.h"
#include "InfiniteHomogenousReactor.h"
#include "HomogenousMesh.h"

HomogenousMicroCell::HomogenousMicroCell(InfiniteHomogenousReactor* reactor, const Real &initial_temperature) : MicroCell(reactor) 
{
    this->_reactor = reactor;
    
    //Solver Settings
    _time_step = _reactor->_input_file_reader->getInputFileParameter("Thermal Time Iteration", static_cast<Real>(100e-9) );
    _start_time = 0;
    _outward_integrated_power = 0;
    _integrated_power = 0;
    _current_time = _start_time;    
    
    //Mesh Setup
    _number_mesh_nodes = 1;
    this->setMesh(new HomogenousMesh(_reactor->_micro_sphere_geometry) );
    _solution.resize( 1 , initial_temperature );
    
    //Boundary Condition   
    _outer_boundary_condition = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();  
   
}

HomogenousMicroCell::~HomogenousMicroCell() {}

void HomogenousMicroCell::setMesh(HomogenousMesh* mesh)
{
    _mesh = mesh;
    MicroCell::setMesh(mesh);
}

void HomogenousMicroCell::getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const
{
    cell_temperature = _solution[0];
    cell_volume = this->_mesh->_volume[0];
}

void HomogenousMicroCell::getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const
{
    cell_temperature = _solution[0];
    cell_volume = this->_mesh->_volume[0];
}

Real HomogenousMicroCell::getVolume() const
{
    Real boundary_volume = 1; // m^3 
    return boundary_volume;
}

MicroSolution HomogenousMicroCell::solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
{
    
    std::vector<Real> previous_solution = this->_solution;
    std::vector<Real> current_solution;    
    current_solution.resize( previous_solution.size() );
    
    Real time = _current_time;
    unsigned long iteration = 0;
    
    while( time < _current_time + simulation_time_step )
    {
       
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        Real dSolution = 0;
        Real outward_heat_flux = 0;
        Real net_heat_flux = 0;
        
        Real temperature = previous_solution[0];
               
        MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getHomogenizedMaterialProperties(temperature);
              
        // Heating Power
        Real internal_power = power_distribution[0];
        _integrated_power += internal_power * _time_step;
            
        //use the boundary condition to determine the outward heat flux
        outward_heat_flux = this->_outer_boundary_condition->getHeatFlux(this, internal_power);
        _outward_integrated_power += -outward_heat_flux * _time_step;
        _outward_current_power = -outward_heat_flux;
            
            
        Real base_thermal_inertia = material_data._density * material_data._specific_heat;
        Real thermal_inertia = base_thermal_inertia;
            
        net_heat_flux = outward_heat_flux + internal_power;  
        dSolution = net_heat_flux / thermal_inertia;
        current_solution[0] = previous_solution[0] + _time_step * dSolution;

        previous_solution = current_solution;
        ++iteration;
        time = time + _time_step;
    }
    
    _current_time = time;
    _solution = current_solution;
    return this->getCurrentMicrosolution();
}

std::vector<Real> HomogenousMicroCell::getRespresentativePowerDistribution(const Real &average_power_density)
{
    std::vector<Real> power_distribution = std::vector<Real>();
    power_distribution.push_back( average_power_density );
    return power_distribution;
}


