/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CompositeMicroCell.cpp
 * Author: chris
 * 
 * Created on March 13, 2017, 1:25 PM
 */
#include <iostream>
#include "CompositeMicroCell.h"

using namespace MaterialLibrary;


CompositeMicroCell::CompositeMicroCell(InfiniteCompositeReactor* reactor, const Real &initial_temperature) : MicroCell(reactor)
{
    this->_reactor = reactor;    
    
    //Solver Settings
    _time_step = _reactor->_input_file_reader->getInputFileParameter("Thermal Time Iteration", static_cast<Real>(100e-9) );
    _start_time = 0;
    _outward_integrated_power = 0;
    _integrated_power = 0;
    _current_time = _start_time;    
    
    //Mesh Setup
    _number_mesh_nodes = _reactor->_input_file_reader->getInputFileParameter("Thermal Mesh Size", static_cast<int>(100) );
    Real min_nodes_per_zone = _reactor->_input_file_reader->getInputFileParameter("Minimum Elements Per Zone", static_cast<int>(6) );
    
    int cells_per_zone = this->_reactor->_input_file_reader->getInputFileParameter("Cells Per Zone", 1 );
    this->setMesh(new SphericalMesh(_reactor->_micro_sphere_geometry,min_nodes_per_zone, _number_mesh_nodes, cells_per_zone));
    _solution.resize( _mesh->numberOfNodes(), initial_temperature );
    
    //Boundary Condition   
    _outer_boundary_condition = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();     
    _solver_order = SolverOrder::SECOND;
}


CompositeMicroCell::~CompositeMicroCell() {}

MicroSolution CompositeMicroCell::solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
{
    if( _solver_order == SolverOrder::SECOND )
    {
        return solveSecondOrder(simulation_time_step,power_distribution);
    }
    else if( _solver_order == SolverOrder::FOURTH )
    {
        return solveFourthOrder(simulation_time_step,power_distribution);
    }    
}

MicroSolution CompositeMicroCell::solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
{
    std::vector<Real> radial_mesh = _mesh->_position;
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
        Real inward_heat_flux = 0; 
        Real net_heat_flux = 0;
        
        for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
        {
            Real temperature = previous_solution[radial_index];                        
            Materials current_material = _mesh->getMaterial(radial_index);
            MaterialDataPacket material_data = MaterialLibrary::getMaterialProperties(current_material,temperature);
            Dimension center_radius = _mesh->getNodeLocation(radial_index);         
            Real internal_power = power_distribution[radial_index] * _mesh->_volume[radial_index];
            _integrated_power += internal_power * _time_step;
            
            
            //Check for boundary conditions
            //If we are at the last node
            if( radial_index == _mesh->numberOfNodes() - 1)
            {
                //use the boundary condition to determine the outward heat flux
                outward_heat_flux = this->_outer_boundary_condition->getHeatFlux(this,inward_heat_flux + internal_power);
                _outward_integrated_power += -outward_heat_flux * _time_step;
                _outward_current_power = -outward_heat_flux;
            }            
            else
            {   
                //Neighbor dependent Properties
                Real outward_neighbor_cell_temperature = previous_solution[radial_index + 1];
                Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
                Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
                Real outward_dTdr = outward_delta_T / outward_distance; 
                
                
                outward_heat_flux = - _mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity;                               
            }   
            
            Real base_thermal_inertia = material_data._density * material_data._specific_heat;
            Real thermal_inertia = base_thermal_inertia * _mesh->_volume[radial_index];
            
            net_heat_flux = inward_heat_flux + outward_heat_flux + internal_power;  
            dSolution = net_heat_flux / thermal_inertia;
            current_solution[ radial_index ] = previous_solution[ radial_index ] + _time_step * dSolution;
            
            //for the next cell the inward heat flux is just the negative of the outward for the previous cell
            inward_heat_flux = -outward_heat_flux;
            
        }        
        
        previous_solution = current_solution;
        ++iteration;
        time = time + _time_step;
    }
    
    _current_time = time;
    _solution = current_solution;
    return this->getCurrentMicrosolution();
}

MicroSolution CompositeMicroCell::solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
{
    std::vector<Real> previous_solution = _solution;
    std::vector<Real> current_solution;
    
    current_solution.resize( previous_solution.size() );
    
    Real time = _current_time;
    long iteration = 0;
    
    while( time < _current_time + simulation_time_step )
    {
        iteration++;
        time = time + _time_step;
        
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        
        
        for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
        {
            Real temperature = previous_solution[radial_index];
            Dimension radial_position = _mesh->getNodeLocation(radial_index); 
            Real power = power_distribution[radial_index];
            MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getMaterialProperties(radial_position,temperature);//testMaterialProperties(radial_position);//
            
            
            /*
            Real dSolution = 0;
            
            //Check for boundary conditions
            if ( radial_index == 0 )
            {    
                 Real d2Tdr2 = (-30 * previous_solution[radial_index] + 32*previous_solution[radial_index+1] - 2 * previous_solution[radial_index + 2])/(12 * _solver_settings._element_size*_solver_settings._element_size );
                dSolution =  ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat );
            } 
            else if( radial_index == _solver_settings._number_points - 1)
            {
                Real d2Tdr2 = ( 32*previous_solution[radial_index - 1] -30*previous_solution[radial_index] - 2 * previous_solution[radial_index ])/(_solver_settings._element_size * _solver_settings._element_size);
                dSolution = ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat );                
            }  
            //if not at boundary calculate as normal
            else if(radial_index == 1)
            {
                Real dTdr = (previous_solution[radial_index+1] - previous_solution[radial_index-1])/(2*_solver_settings._element_size);
                Real d2Tdr2 = ( previous_solution[radial_index+1] - 2*previous_solution[radial_index] + previous_solution[radial_index-1] )/(_solver_settings._element_size*_solver_settings._element_size);
                dSolution =  ( material_data._thermal_conductivity_temperature_derivative * dTdr * dTdr * radial_position*radial_position + material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2) + power*radial_position*radial_position)/(radial_position*radial_position*( material_data._density * material_data._specific_heat + material_data._specific_heat_temperature_derivative * material_data._density ) );
            }
            else if( radial_index == _solver_settings._number_points - 2 )
            {
                Real dTdr = (previous_solution[radial_index+1] - previous_solution[radial_index-1])/(2*_solver_settings._element_size);
                Real d2Tdr2 = ( previous_solution[radial_index+1] - 2*previous_solution[radial_index] + previous_solution[radial_index-1] )/(_solver_settings._element_size*_solver_settings._element_size);
                dSolution =  ( material_data._thermal_conductivity_temperature_derivative * dTdr * dTdr * radial_position*radial_position + material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2) + power*radial_position*radial_position)/(radial_position*radial_position*( material_data._density * material_data._specific_heat + material_data._specific_heat_temperature_derivative * material_data._density ) );               
            }
            else
            {
                Real dTdr = (-previous_solution[radial_index+2] + 8*previous_solution[radial_index+1]  - 8*previous_solution[radial_index-1] + previous_solution[radial_index-2] )/(12*_solver_settings._element_size);
                //Real d2Tdr2 = ( previous_solution[radial_index+1] - 2*previous_solution[radial_index] + previous_solution[radial_index-1] )/(_solver_settings._element_size*_solver_settings._element_size);
                
                Real d2Tdr2 = ( -1 * previous_solution[radial_index+2] + 16*previous_solution[radial_index + 1] - 30*previous_solution[radial_index] + 16*previous_solution[radial_index-1] - previous_solution[radial_index - 2] )/(12 * _solver_settings._element_size*_solver_settings._element_size);
                //std::cout<<d2Tdr2 << " " << d2Tdr22 << "\n";
                
                dSolution =  ( material_data._thermal_conductivity_temperature_derivative * dTdr * dTdr * radial_position*radial_position + material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2) + power*radial_position*radial_position)/(radial_position*radial_position*( material_data._density * material_data._specific_heat + material_data._specific_heat_temperature_derivative * material_data._density ) ); 
            }
            
            current_solution[ radial_index ] = previous_solution[ radial_index ] + _solver_settings._time_step * dSolution;
            */
            
        }        
        
        previous_solution = current_solution;
    }
    
    this->_current_time = time;
    this->_solution = current_solution;    
    return this->getCurrentMicrosolution();
}


std::vector<MicroSolution> CompositeMicroCell::iterateInitialConditions(const std::vector<Real> &power_distribution)
{
    MicroSolution initial_solution = this->presolveSteadyStateAnalytical(power_distribution);
    
    std::vector<MicroSolution> solutions = std::vector<MicroSolution>();
    solutions.push_back(initial_solution);
    
    Real solution_time_step = _time_step * 50000;
    Real max_residual = 1;
    Real avg_residual = 0;
    Real index = 0;
    Real desired_residual = _reactor->_input_file_reader->getInputFileParameter("Steady State Temperature Solution Max Residual", static_cast<Real>(0.005) );
    
    
    while( max_residual > desired_residual )
    {
       MicroSolution solution = this->solve(solution_time_step,power_distribution); 
       
       if(index > 0)
       {
            //calculate the maximum change in solution from previous to new 
            vector_residuals(solution._solution, solutions.back()._solution, max_residual, avg_residual);
            std::cout<<"Residual: " << max_residual << " Time Step: " << solution_time_step << std::endl;
       }
       solutions.push_back(solution);
       
       index++;
       
    }
       
    return solutions;
}


MicroSolution CompositeMicroCell::presolveSteadyStateAnalytical(const std::vector<Real> &power_distribution)
{
    for( int radial_index = _number_mesh_nodes - 2; radial_index >= 0; --radial_index )
    {
        
        //Averaged material properties in the cell assuming right side temperature material properties
        Real temperature_to_the_right =  _solution[radial_index + 1];
        MaterialDataPacket material_data = MaterialLibrary::getMaterialProperties(_mesh->getMaterial(radial_index) ,temperature_to_the_right);
        
        Real total_power = 0;
        
        for( int power_index = 0; power_index <= radial_index; power_index++ )
        {
             Real volume = _mesh->_volume[power_index];
             total_power += volume * power_distribution[power_index];
        }
       
        Real cell_outer_surface_area = _mesh->_outer_surface[radial_index];
        Real outer_surface_heat_flux = total_power / cell_outer_surface_area;        
        Real power_density = power_distribution[radial_index];    
        Dimension right_radial_position = _mesh->_outer_radius[radial_index];
        Dimension left_radial_position = _mesh->_inner_radius[radial_index];
        
        const Real constant_1 = right_radial_position*right_radial_position * (outer_surface_heat_flux / material_data._thermal_conductivity - power_density * right_radial_position / (3.0 * material_data._thermal_conductivity ) );
        const Real constant_2 = temperature_to_the_right + power_density * right_radial_position * right_radial_position / ( 6.0 * material_data._thermal_conductivity) - constant_1 / right_radial_position;

        if( left_radial_position > 0)
        {
            _solution[radial_index] = -power_density * left_radial_position*left_radial_position / (6.0 * material_data._thermal_conductivity) + constant_1/left_radial_position + constant_2; 
        }
        else
        {
             _solution[radial_index] = -power_density * left_radial_position*left_radial_position / (6.0 * material_data._thermal_conductivity) + constant_2; 
        }

    }      
    return MicroSolution(_mesh->_position, _solution,0);
       
}

void CompositeMicroCell::setOuterMaterialTemperature(const Real &outer_temperature)
{
    for(int index = 0; index < _mesh->numberOfNodes(); ++index)
    {
        if(_mesh->_zone[index] != 1)
        {
            _solution[index] = outer_temperature;
        }
    }
}

std::vector<Real> CompositeMicroCell::getRespresentativePowerDistribution(const Real &average_power_density)
{
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    Dimension number_points = _mesh->numberOfNodes();
    power_distribution.reserve( number_points);    
    
    //Get the boundary locations     
    Real fuel_kernel_radius = _reactor->_micro_sphere_geometry->getFuelKernelRadius();
    Real half_box_edge = _reactor->_micro_sphere_geometry->getOuterRadius();
    Real kernel_power = average_power_density * 8 * pow(half_box_edge,3)/sphere_volume(fuel_kernel_radius);
    Real power = 0;
    
    
    for( size_t index = 0; index < number_points; ++index)
    {
        if( _mesh->_zone[index] == 1  )
        {
            power = kernel_power;           
        }
        else
        {
            power = 0;           
        }
        
        power_distribution.push_back( power );
    }
    
    return power_distribution;
    
    
    
}

void CompositeMicroCell::getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const
{
    Dimension material_inner_radius = 0.0;
    Dimension material_outer_radius = _reactor->_micro_sphere_geometry->_geometry[zone].second;
    
    if( zone > 0 )
    {
        material_inner_radius = _reactor->_micro_sphere_geometry->_geometry[zone-1].second;
    }
    
    int number_zones = _reactor->_micro_sphere_geometry->_geometry.size();  
    
    Dimension outer_radius;
    Dimension material_delta_radius = material_outer_radius - material_inner_radius;
    Dimension inner_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division - 1)/static_cast<Real>(zone_divisions) );
        
    
    //If we are in the last zone and last zone division we make the radius the equivalent to the neutronics square volume
    if( zone == number_zones -1 && zone_divisions == current_division)
    {
        Dimension a = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division )/static_cast<Real>(zone_divisions) );
        
        outer_radius =a * std::pow( 6.0/M_PI , 1.0/3.0);
    }
    else
    {
        outer_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division )/static_cast<Real>(zone_divisions) );
    }
    this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature, cell_volume);
}

/**
 * Check to see if this mesh at this mesh index is inside the outer and inner radius
 * @param mesh_index
 * @param inner_radius
 * @param outer_radius
 * @return 
 */
inline bool CompositeMicroCell::inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
{
    return  _mesh->_inner_radius[mesh_index] < outer_radius && _mesh->_outer_radius[mesh_index] > inner_radius;
}

inline Real CompositeMicroCell::overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
{
    if( inMeshCell( mesh_index, inner_radius, outer_radius) )
    {
        //if it is completely in the cell
        if( _mesh->_inner_radius[mesh_index] < inner_radius && _mesh->_outer_radius[mesh_index] > outer_radius )
        {
            return _mesh->_volume[mesh_index];
        }
        //The cell over engulfs the cell
        else if( _mesh->_inner_radius[mesh_index] > inner_radius && _mesh->_outer_radius[mesh_index] < outer_radius )
        {
            return _mesh->_volume[mesh_index];
        }
        //There is overlap between regions
        else 
        {
            Real inner_r, outer_r;
            
            //There is left side overlap
            if( _mesh->_inner_radius[mesh_index] < inner_radius )
            {
                inner_r = inner_radius;
            }
            else
            {
                inner_r = _mesh->_inner_radius[mesh_index];
            }
            
            //There is right side overlap
            if(  _mesh->_outer_radius[mesh_index] > outer_radius )
            {
                outer_r = outer_radius;
            }
            else
            {
                outer_r = _mesh->_outer_radius[mesh_index];
            }
            
            return (sphere_volume(outer_r) - sphere_volume(inner_r));
        }
    }
    else
    {
        return 0;
    }
}

void CompositeMicroCell::getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume) const
{
    cell_volume = sphere_volume(outer_radius) - sphere_volume(inner_radius);    
    Real integrated_volume = 0;
    Real temperature_weighted_volume = 0;
    
    for(int index = 0; index < _mesh->numberOfNodes(); index++  )
    {
        Real unit_mesh_volume =  overlappingVolumeWithMeshCell( index, inner_radius, outer_radius);
        temperature_weighted_volume += _solution[index] * unit_mesh_volume;
        integrated_volume += unit_mesh_volume;
    }
    
    cell_temperature = temperature_weighted_volume/cell_volume;
    
}

void CompositeMicroCell::getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const
{
    Dimension inner_radius = 0.0;
    Dimension outer_radius = _reactor->_micro_sphere_geometry->_geometry[zone].second;
    
    if( zone > 0 )
    {
        inner_radius = _reactor->_micro_sphere_geometry->_geometry[zone-1].second;
    }
    
   return this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature,cell_volume);
}

void CompositeMicroCell::setMesh(SphericalMesh* mesh)
{
    this->_mesh = mesh;
    this->MicroCell::setMesh(mesh);
}