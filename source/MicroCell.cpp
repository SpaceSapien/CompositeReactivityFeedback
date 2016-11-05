/*
 * The MicroCell is a spherical particle
 * 
 */

/* 
 * File:   MicroHeatTransferExplicitSolver.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 5:46 PM
 */
//#define VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL
#include <iostream>
#include <vector>
#include <cmath>
#include "MicroSolution.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "MaterialDataPacket.h"
#include "InputDataFunctions.h"
#include "MicroCellBoundaryCondition.h"
#include "MicroCell.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"
#include "RadialMesh.h"



MicroCell::MicroCell() {}

MicroCell::MicroCell(InfiniteCompositeReactor* reactor,const Real &initial_temperature) 
{
    
    this->_reactor = reactor;    
    
    //Solver Settings
    _time_step = _reactor->_input_file_reader->getInputFileParameter("Thermal Time Iteration", static_cast<Real>(100e-9) );
    _solver_order = SolverOrder::SECOND;
    _start_time = 0;
    _current_time = _start_time;    
    
    //Mesh Setup
    _number_mesh_nodes = _reactor->_input_file_reader->getInputFileParameter("Thermal Mesh Size", static_cast<int>(100) );
    Real min_nodes_per_zone = _reactor->_input_file_reader->getInputFileParameter("Minimum Elements Per Zone", static_cast<int>(6) );
    _mesh = new RadialMesh(_reactor->_micro_sphere_geometry,min_nodes_per_zone, _number_mesh_nodes);
    _solution.resize( _mesh->numberOfNodes(), initial_temperature );
    
    //Boundary Condition   
    _outer_boundary_condition = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();    
    
}

MicroCell::~MicroCell()
{
    delete _mesh;
    delete _outer_boundary_condition;
}


MicroSolution MicroCell::presolveSteadyStateAnalytical(const Real &average_power_density)
{

    std::vector<Real> power_distribution = this->getRepresentativeKernelPowerDistribution(average_power_density);
    
    for( int radial_index = _number_mesh_nodes - 2; radial_index >= 0; --radial_index )
    {
        Dimension radial_position = _mesh->getNodeLocation(radial_index);
        
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

MicroSolution MicroCell::solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
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

MicroSolution MicroCell::solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
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


MicroSolution MicroCell::getCurrentMicrosolution()
{
    MicroSolution solution = MicroSolution(_mesh->_position,_solution,_current_time);
    return solution;
}

void MicroCell::setBoundaryCondition(MicroCellBoundaryCondition* boundary_condition)
{
    if(boundary_condition != nullptr)
    {
        delete _outer_boundary_condition;
    }
    
    this->_outer_boundary_condition = boundary_condition;
}

MicroSolution MicroCell::solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
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
            
            //Check for boundary conditions
            //If we are at the last node
            if( radial_index == _mesh->numberOfNodes() - 1)
            {
                //no outward heat flux
                outward_heat_flux = this->_outer_boundary_condition->getHeatFlux(this,inward_heat_flux + internal_power);
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
            
            Real thermal_inertia = (material_data._density * material_data._specific_heat  + temperature * ( material_data._density * material_data._specific_heat_temperature_derivative + material_data._specific_heat * material_data._density_temperature_derivative) ) * _mesh->_volume[radial_index];
            
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


MaterialDataPacket MicroCell::testMaterialProperties(const Real &radius)
{
    Real density, specific_heat, thermal_conductivity;
    
    if(radius < 2e-4 )
    {   
        density = 1.0928e+04;
        specific_heat = 298.4734;
        thermal_conductivity = 6.9638;
    } 
    else
    {
        density = 1.7069e+03;
        specific_heat = 982;
        thermal_conductivity = 93.5000;
    }   

    MaterialDataPacket mdp = MaterialDataPacket(thermal_conductivity, density, specific_heat,0,0,0);
    
    return mdp;
    
}


std::vector<MicroSolution> MicroCell::iterateInitialConditions(const Real &initial_average_power_density)
{
    MicroSolution initial_solution = this->presolveSteadyStateAnalytical(initial_average_power_density);
    
    std::vector<Real> power_distribution =  this->getRespresentativePowerDistribution(initial_average_power_density);    
    std::vector<MicroSolution> solutions = std::vector<MicroSolution>();
    solutions.push_back(initial_solution);
    
    Real solution_time_step = 0.05;
    Real max_residual = 1;
    Real index = 0;
    Real desired_residual = _reactor->_input_file_reader->getInputFileParameter("Steady State Temperature Solution Max Residual", static_cast<Real>(0.005) );
    
    
    while( max_residual > desired_residual )
    {
       MicroSolution solution = MicroCell::solve(solution_time_step,power_distribution); 
       
       if(index > 0)
       {
            //calculate the maximum change in solution from previous to new 
            max_residual =  calculationMaximumResidual(solution._solution, solutions.back()._solution, solution_time_step);
            std::cout<<"Residual: " << max_residual << " Time Step: " << solution_time_step << std::endl;
       }
       solutions.push_back(solution);
       
       index++;
       
    }
       
    return solutions;
}

Real MicroCell::calculationMaximumResidual(const std::vector<Real> &vector_1, const std::vector<Real> &vector_2, const Real &time_step)
{
    size_t vector_1_size = vector_1.size();
    size_t vector_2_size = vector_2.size();
    
    if( vector_1_size == vector_2_size )
    {
        Real maximum_residual = 0;    
        
        for(int index = 0; index < vector_1_size; index++)
        {
            Real index_residual = std::abs( 1.0 - 1.0 * vector_2[index]/vector_1[index] );
            
            if( index_residual > maximum_residual)
            {
                maximum_residual = index_residual;
            }
        }
        
        return maximum_residual;
    }
    else
    {
        return -1;
    }
}

void MicroCell::getCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume )
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
 
inline bool MicroCell::inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
{
    return  _mesh->_inner_radius[mesh_index] < outer_radius && _mesh->_outer_radius[mesh_index] > inner_radius;
}

inline Real MicroCell::overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
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

void MicroCell::getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume)
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

void MicroCell::getAverageTemperature(const int &zone, Real &cell_temperature, Real &cell_volume)
{
    Dimension inner_radius = 0.0;
    Dimension outer_radius = _reactor->_micro_sphere_geometry->_geometry[zone].second;
    
    if( zone > 0 )
    {
        inner_radius = _reactor->_micro_sphere_geometry->_geometry[zone-1].second;
    }
    
   return this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature,cell_volume);
}

std::vector<Real> MicroCell::getRespresentativePowerDistribution(const Real &average_power_density)
{
    #ifdef HOMOGENIZE_POWER

    return this->getRespresentativeHomogenizedPowerDistribution(average_power_density);
    
    #else

    return this->getRepresentativeKernelPowerDistribution(average_power_density);
    
    #endif
}

std::vector<Real> MicroCell::getRepresentativeKernelPowerDistribution(const Real &average_power_density)
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
        if( _mesh->_zone[index] == 0  )
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


std::vector<Real> MicroCell::getRespresentativeHomogenizedPowerDistribution(const Real &average_power_density)
{
    //Get the radial mesh
    size_t number_points = _mesh->numberOfNodes();
    
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    power_distribution.reserve( number_points);
    
    for( size_t index = 0; index < number_points; ++index)
    {
        power_distribution.push_back( average_power_density );
    }
    
    return power_distribution;
}
