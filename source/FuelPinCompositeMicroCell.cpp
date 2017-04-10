/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinCompositeMicroCell.cpp
 * Author: chris
 * 
 * Created on April 6, 2017, 6:18 PM
 */

#include "FuelPinCompositeMicroCell.h"
#include "ReactorKinetics.h"

FuelPinCompositeMicroCell::FuelPinCompositeMicroCell(FuelPinReactor* reactor, const Real &temperature,const Real &macro_scale_position) : CompositeMicroCell(reactor, temperature)
{
    _macro_scale_position = macro_scale_position;
   // _time_step /= 10.0;
}
FuelPinCompositeMicroCell::~FuelPinCompositeMicroCell() {}

std::vector<MicroSolution> FuelPinCompositeMicroCell::iterateInitialConditions(const std::vector<Real> &power_distribution)
{
    std::vector<MicroSolution> return_solutions = this->CompositeMicroCell::iterateInitialConditions(power_distribution);
    
    Real initial_power = _mesh->getVolumeWeightedQuantity(power_distribution);
    
    this->solveAndSetSurfaceResistance(initial_power);
    return return_solutions;
}
void FuelPinCompositeMicroCell::solveAndSetSurfaceResistance(const Real &initial_power_per_cell)
{
    Real matrix_avg_temperature, matrix_volume;
    this->getAverageZoneTemperature(this->_number_zones, matrix_avg_temperature, matrix_volume );
    
    Real coating_radius = _reactor->_micro_sphere_geometry->_geometry[_number_zones-2].second;
    _coating_surface_area = 4.0 * M_PI *  coating_radius * coating_radius;
    Real heat_per_unit_area = initial_power_per_cell / _coating_surface_area; 
    
    
    int number_materials = _reactor->_micro_sphere_geometry->_geometry.size();
    Real outer_coating_radius = _reactor->_micro_sphere_geometry->_geometry[number_materials - 2].second;
    Real matrix_box_half_length = _reactor->_micro_sphere_geometry->_geometry[number_materials - 1].second;    
    Real total_volume = 8.0 * matrix_box_half_length*matrix_box_half_length*matrix_box_half_length;
    Real inner_volume = sphere_volume(outer_coating_radius);
    _number_particles_per_volume = 1.0 / total_volume;
    _matrix_volume_fraction = (total_volume - inner_volume) / total_volume;
    
    for( std::size_t node_index = 0; node_index < _mesh->numberOfNodes(); ++node_index)
    {
        if(_mesh->_zone[ node_index ] == _number_zones)
        {
            Real coating_temperature = _solution[node_index];
            Real delta_t = coating_temperature - matrix_avg_temperature;
            
            _surface_resistance = heat_per_unit_area / delta_t;
            return;               
        }
    }
    
    
   
}


Real FuelPinCompositeMicroCell::solveMicroCellCoupling(const Real &current_macrocell_temperature, const Real &microcell_internal_power_density, const Real &time_step)
{
    std::vector<Real> current_solution;
    current_solution.resize( _solution.size() );
    
    Real micro_scale_outward_power;
    //Smear the power over the center of the cell
    std::vector<Real> power_data = this->getRespresentativePowerDistribution(microcell_internal_power_density);
    bool break_loop = false;
    
    Real outward_heat_flux = 0;
    Real inward_heat_flux = 0; 

    for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
    {
        
        
        Real temperature = _solution[radial_index];                        
        Materials current_material = _mesh->getMaterial(radial_index);
        MaterialDataPacket material_data = MaterialLibrary::getMaterialProperties(current_material,temperature);
        Dimension center_radius = _mesh->getNodeLocation(radial_index);         
        Real internal_power = power_data[radial_index] * _mesh->_volume[radial_index];
        _integrated_power += internal_power * time_step;


        //Check for boundary conditions
        //If we are at the boundary between then matrix and the outermost coating/kernel
        if( _mesh->_zone[radial_index] == this->_number_zones)
        {
            Real delta_t = (_solution[radial_index] - current_macrocell_temperature);
            //use the boundary condition to determine the outward heat flux
            outward_heat_flux = -_coating_surface_area * _surface_resistance * delta_t;  
            //number_particles_per_volume *matrix_volume_fraction) * ;
            Real volumetric_power = _number_particles_per_volume * outward_heat_flux;
            micro_scale_outward_power = volumetric_power;
            _outward_integrated_power += -outward_heat_flux * time_step;
            _outward_current_power = -outward_heat_flux;
            break_loop = true;
        }
        else
        {   
            //Neighbor dependent Properties
            Real outward_neighbor_cell_temperature = _solution[radial_index + 1];
            Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
            Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
            Real outward_dTdr = outward_delta_T / outward_distance; 
            outward_heat_flux = - _mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity;                               
        }   

        Real base_thermal_inertia = material_data._density * material_data._specific_heat;
        Real thermal_inertia = base_thermal_inertia * _mesh->_volume[radial_index];

        Real net_heat_flux = inward_heat_flux + outward_heat_flux + internal_power;  
        Real dSolution = net_heat_flux / thermal_inertia;
        current_solution[ radial_index ] = _solution[ radial_index ] + time_step * dSolution;

        //for the next cell the inward heat flux is just the negative of the outward for the previous cell
        inward_heat_flux = -outward_heat_flux;

        if( break_loop )
        {
            for( int second_radial_index = radial_index + 1; second_radial_index < _mesh->numberOfNodes(); ++second_radial_index)
            {
                current_solution[ second_radial_index ] = current_macrocell_temperature;
            }
            break;
        }

    }        

           
    _current_time = _current_time + time_step;
    _solution = current_solution;
    return micro_scale_outward_power;
}