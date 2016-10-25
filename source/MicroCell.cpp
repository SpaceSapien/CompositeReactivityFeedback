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
#define VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL
#include <iostream>
#include <vector>
#include <cmath>
#include "MicroSolution.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "ExplicitSolverSettings.h"
#include "MaterialDataPacket.h"
#include "InputDataFunctions.h"
#include "MicroCellBoundaryCondition.h"
#include "MicroCell.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"

MicroCell::MicroCell() {}

MicroCell::MicroCell(InfiniteCompositeReactor* reactor,const Real &initial_temperature) 
{
    
    this->_reactor = reactor;
    Dimension max_radius = _reactor->_micro_sphere_geometry->getOuterRadius();
    this->_solver_settings = ExplicitSolverSettings(ExplicitSolverSettings::SolverOrder::SECOND, max_radius );     
    
    
    //Default initialization to 350 K
    this->_solution.resize( _solver_settings._radial_mesh.size(), initial_temperature );
    
    this->_current_time = _solver_settings._start_time;
    this->_outer_boundary_condition = MicroCellBoundaryCondition::getRefelectedBoundaryCondition();
}

MicroSolution MicroCell::solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
{
    if( _solver_settings._solver_order == ExplicitSolverSettings::SolverOrder::SECOND )
    {
        return solveSecondOrder(simulation_time_step,power_distribution);
    }
    else if( _solver_settings._solver_order == ExplicitSolverSettings::SolverOrder::FOURTH )
    {
        return solveFourthOrder(simulation_time_step,power_distribution);
    }    
}

MicroSolution MicroCell::solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
{
    std::vector<Real> radial_mesh = _solver_settings._radial_mesh;
    std::vector<Real> previous_solution = _solution;
    std::vector<Real> current_solution;
    
    current_solution.resize( previous_solution.size() );
    
    Real time = _current_time;
    long iteration = 0;
    
    while( time < _current_time + simulation_time_step )
    {
        iteration++;
        time = time + _solver_settings._time_step;
        
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        
        
        for( int radial_index = 0; radial_index < _solver_settings._number_points; ++radial_index )
        {
            Real temperature = previous_solution[radial_index];
            Dimension radial_position = radial_mesh[radial_index]; 
            Real power = power_distribution[radial_index];
            MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getMaterialProperties(radial_position,temperature);//testMaterialProperties(radial_position);//
            
            
            
            Real dSolution = 0;
            
            //Check for boundary conditions
            if ( radial_index == 0 )
            {    
                 Real d2Tdr2 = (-30*previous_solution[radial_index] + 32*previous_solution[radial_index+1] - 2 * previous_solution[radial_index + 2])/(12 * _solver_settings._element_size*_solver_settings._element_size );
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
            
            
        }        
        
        previous_solution = current_solution;
    }
    
    this->_current_time = time;
    this->_solution = current_solution;    
    return this->getCurrentMicrosolution();
}


MicroSolution MicroCell::getCurrentMicrosolution()
{
    MicroSolution solution = MicroSolution(_solver_settings._radial_mesh,_solution,_current_time);
    return solution;
}

void MicroCell::setBoundaryCondition(const MicroCellBoundaryCondition &boundary_condition)
{
    std::vector<Real> radial_mesh = _solver_settings._radial_mesh;
    this->_outer_boundary_condition = boundary_condition;
}

MicroSolution MicroCell::solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution)
{
    
    std::vector<Real> radial_mesh = _solver_settings._radial_mesh;
    std::vector<Real> previous_solution = this->_solution;
    std::vector<Real> current_solution;
    
    current_solution.resize( previous_solution.size() );
    
    Real time = _current_time;
    long iteration = 0;
    
    while( time < _current_time + simulation_time_step )
    {
        iteration++;
        time = time + _solver_settings._time_step;
        
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        for( int radial_index = 0; radial_index < _solver_settings._number_points; ++radial_index )
        {
            Real temperature = previous_solution[radial_index];
            Dimension radial_position = radial_mesh[radial_index]; 
            Real power = power_distribution[radial_index];
            MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getMaterialProperties(radial_position,temperature); //testMaterialProperties(radial_position);
            Real dSolution = 0;
            Real dTdr, d2Tdr2;
            
            //Check for boundary conditions
            if ( radial_index == 0 )
            {    
                d2Tdr2 = ( 2*previous_solution[radial_index + 1] - 2*previous_solution[radial_index] )/(_solver_settings._element_size*_solver_settings._element_size);
                
                #ifdef VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL

                    dSolution =  ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat + material_data._density * material_data._thermal_conductivity_temperature_derivative * temperature );
                    
                #else
                    
                    dSolution =  ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat );
                    
                #endif
            }   
            else if( radial_index == _solver_settings._number_points - 1)
            {
                dTdr = _outer_boundary_condition.getdTdr(_solver_settings._solver_order, -1,  previous_solution , _solver_settings._radial_mesh);
                d2Tdr2 = _outer_boundary_condition.getd2Tdr2(_solver_settings._solver_order, -1,  previous_solution, _solver_settings._radial_mesh);
                
                
                #ifdef VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL

                    dSolution =  ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat + material_data._density * material_data._thermal_conductivity_temperature_derivative * temperature );
                    
                #else
                    
                    dSolution =  ( material_data._thermal_conductivity *(d2Tdr2) + power)/( material_data._density * material_data._specific_heat );
                    
                #endif
            }  
            //if not at boundary calculate as normal
            else
            {
                Real dTdr = (previous_solution[radial_index+1] - previous_solution[radial_index-1])/(2*_solver_settings._element_size);
                Real d2Tdr2 = ( previous_solution[radial_index+1] - 2*previous_solution[radial_index] + previous_solution[radial_index-1] )/(_solver_settings._element_size*_solver_settings._element_size);
                
                #ifdef VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL

                    dSolution =  ( material_data._thermal_conductivity_temperature_derivative*dTdr*dTdr*radial_position*radial_position + material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2) + power*radial_position*radial_position)/(radial_position*radial_position * ( material_data._density * material_data._specific_heat + material_data._density * material_data._specific_heat_temperature_derivative * temperature ) );
                    /*std::cout<<radial_position<<"\n";
                    std::cout<<material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2)<<"\n";
                    std::cout<<material_data._thermal_conductivity_temperature_derivative*dTdr*dTdr*radial_position*radial_position << "\n";
                    std::cout<<power*radial_position*radial_position << "\n";
                    std::cout<<material_data._density * material_data._specific_heat << "\n";
                    std::cout<<material_data._density * material_data._specific_heat_temperature_derivative * temperature << "\n";
                    std::cout<<"\n";*/
                    
                #else
                    
                    dSolution =  (material_data._thermal_conductivity*(2*radial_position*dTdr + radial_position*radial_position*d2Tdr2) + power*radial_position*radial_position)/(radial_position*radial_position * material_data._density * material_data._specific_heat  );
                              
                #endif                
                
            }
            
            current_solution[ radial_index ] = previous_solution[ radial_index ] + _solver_settings._time_step * dSolution;
            
            
        }        
        
        previous_solution = current_solution;
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
    std::vector<Real> power_distribution =  this->getRespresentativePowerDistribution(initial_average_power_density);    
    std::vector<MicroSolution> solutions = std::vector<MicroSolution>();
    Real solution_time_step = 0.05;
    Real max_residual = 1;
    Real index = 0;
    Real desired_residual = _reactor->_input_file_reader->getInputFileParameter("Steady State Temperature Solution Max Residual", static_cast<Real>(0.01) );
    
    
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
    
    Dimension material_delta_radius = material_outer_radius - material_inner_radius;
    Dimension inner_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division - 1)/static_cast<Real>(zone_divisions) );
    Dimension outer_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division )/static_cast<Real>(zone_divisions) );
    
    this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature, cell_volume);
}
 
void MicroCell::getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume)
{
    Real element_size = _solver_settings._element_size;
    Real volume_weighted_temperature = 0;
    Real total_volume = 0;    
    std::vector<Real> radial_mesh = this->_solver_settings._radial_mesh;
    
    for(int index = 0; index < radial_mesh.size(); index++  )
    {
        //If the grid point is fully within the cell
        if( radial_mesh[index] < outer_radius - element_size/2 && radial_mesh[index] - element_size/2 > inner_radius )
        {
            Real volume_slice;

            //The first slice of volume is special             
            if( index == 0)
            {
                volume_slice = 4.0/3.0 * M_PI * pow( element_size/2 , 3);            
            }
            else
            {
                volume_slice = 4.0/3.0 * M_PI * ( pow( radial_mesh[index] + element_size/2 , 3) - pow( radial_mesh[index] - element_size/2 , 3) );
            }

            total_volume += volume_slice;
            Real temperature = _solution[index];
            volume_weighted_temperature += volume_slice * temperature;
        }        
        //If the grid point is only partially in the cell on the outside
        else if( (radial_mesh[index] > outer_radius && radial_mesh[index] - element_size/2 < outer_radius) || (radial_mesh[index] < outer_radius && radial_mesh[index] + element_size/2 > outer_radius) )
        {
            Real volume_slice;
            Real smaller_element_size = outer_radius - ( radial_mesh[index] - element_size/2 );
            
            volume_slice = 4.0/3.0 * M_PI * ( pow( outer_radius  , 3) - pow( radial_mesh[index] - element_size/2 , 3) );
                      
            total_volume += volume_slice;
            Real temperature = _solution[index];
            volume_weighted_temperature += volume_slice * temperature;
        }
        //If the grid point is only partially in the cell on the inside
        else if( ( radial_mesh[index] < inner_radius && radial_mesh[index] + element_size/2 > inner_radius ) || ( radial_mesh[index] > inner_radius && radial_mesh[index] - element_size / 2 < inner_radius ) )
        {
            Real volume_slice;
            Real smaller_element_size = radial_mesh[index] + element_size/2 - inner_radius;
            
            volume_slice = 4.0/3.0 * M_PI * ( pow( radial_mesh[index] + element_size/2 , 3) - pow( inner_radius , 3) );
            
            total_volume += volume_slice;
            Real temperature = _solution[index];
            volume_weighted_temperature += volume_slice * temperature;
        }
    }
    
    cell_temperature = volume_weighted_temperature/total_volume;
    cell_volume = total_volume;
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
    //Get the radial mesh
    std::vector<Real> radial_mesh = _solver_settings._radial_mesh;
    size_t number_points = radial_mesh.size();
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    power_distribution.reserve( number_points);
    //Get the boundary locations     
    Real fuel_kernel_radius = _reactor->_micro_sphere_geometry->getFuelKernelRadius();
    Real outer_radius = _reactor->_micro_sphere_geometry->getOuterRadius();
    
    Real kernel_power = average_power_density * pow(outer_radius,3)/pow(fuel_kernel_radius,3);
    
    for( size_t index = 0; index < number_points; ++index)
    {
        Dimension radial_point = radial_mesh[index];
        Real power;
        
        if(radial_point < fuel_kernel_radius )
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
    std::vector<Real> radial_mesh = _solver_settings._radial_mesh;
    size_t number_points = radial_mesh.size();
    
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    power_distribution.reserve( number_points);
    
    for( size_t index = 0; index < number_points; ++index)
    {
        power_distribution.push_back( average_power_density );
    }
    
    return power_distribution;
}

Real MicroCell::getOuterDerivative()
{
    size_t solution_size_index = _solution.size() -1;
    //use a double spacing limit
    return (_solution[solution_size_index] - _solution[solution_size_index -2])/(2 * _solver_settings._element_size);
}