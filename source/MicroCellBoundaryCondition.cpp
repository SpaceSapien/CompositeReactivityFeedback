/*
 * 
 */

/* 
 * File:   MicroCellBoundaryCondition.cpp
 * Author: chris
 * 
 * Created on December 17, 2015, 4:22 PM
 */

#include "MicroCellBoundaryCondition.h"

MicroCellBoundaryCondition::MicroCellBoundaryCondition() { }


MicroCellBoundaryCondition::MicroCellBoundaryCondition(const BoundaryType &boundary_type) 
{
    _boundary = boundary_type;
}


MicroCellBoundaryCondition  MicroCellBoundaryCondition::getRefelectedBoundaryCondition()
{
    MicroCellBoundaryCondition boundary_condition = MicroCellBoundaryCondition(BoundaryType::Reflected);
    return boundary_condition; 
}

/**
 * Factory function that makes a new boundary conidtion
 * 
 * @param fixed_spatial_derivative dr/dT = cost where const is fixed_spatial_derivative
 * @return MicroCellBoundaryCondition
 */
MicroCellBoundaryCondition  MicroCellBoundaryCondition::getFixedDerivativeBoundaryCondition(const Real &fixed_spatial_derivative)
{
    MicroCellBoundaryCondition boundary_condition = MicroCellBoundaryCondition(BoundaryType::FixedDerivative);
    boundary_condition._fixed_temperature_derivative = fixed_spatial_derivative;
    return boundary_condition; 
}

    
MicroCellBoundaryCondition MicroCellBoundaryCondition::getFixedBoundaryCondition(const Real &fixed_temperature)
{
    MicroCellBoundaryCondition boundary_condition = MicroCellBoundaryCondition(BoundaryType::Fixed);
    boundary_condition._fixed_temperature = fixed_temperature;
    return boundary_condition; 
}

Real MicroCellBoundaryCondition::getdTdr(const ExplicitSolverSettings::SolverOrder &order,const int &index,const std::vector<Real> &previous_solution, const std::vector<Real> &grid)
{
    size_t boundary_index = previous_solution.size() - 1;
    Real delta_r = grid[boundary_index] - grid[boundary_index - 1];
    
    if(order == ExplicitSolverSettings::SolverOrder::SECOND)
    {
        switch(_boundary)
        {
            case BoundaryType::Fixed :
            {
                return  ( _fixed_temperature - previous_solution[boundary_index - 1] ) / ( 2 * delta_r);
                break;
            }
            case BoundaryType::Reflected :
            {
                return 0;
                break;
            }
            case BoundaryType::FixedDerivative :
            {
                return  -1 * _fixed_temperature_derivative ;
                break;
            }
            default :
            {
                throw Error::UnknownBoundaryType;
            }
        }
    }
    //Fourth order spatial derivatives
    else if(order == ExplicitSolverSettings::SolverOrder::FOURTH)
    {
        //one away from right hand boundary
        if(index == -1)
        {
            
        }
        //one away from right hand boundary
        if(index == -2)
        {
            
        }
    }
}

Real MicroCellBoundaryCondition::getd2Tdr2(const ExplicitSolverSettings::SolverOrder &order,const int &index,const std::vector<Real> &previous_solution, const std::vector<Real> &grid )
{
    size_t boundary_index = previous_solution.size() - 1;
    
    Real delta_r = grid[boundary_index] - grid[boundary_index - 1];
    
    if(order == ExplicitSolverSettings::SolverOrder::SECOND)
    {
        
        switch(_boundary)
        {
            case BoundaryType::Fixed :
            {
                return  (previous_solution[boundary_index - 1] - 2 * previous_solution[boundary_index] + _fixed_temperature) / ( delta_r * delta_r );
                break;
            }
            case BoundaryType::Reflected :
            {
                return  (2 * previous_solution[boundary_index - 1] - 2 * previous_solution[boundary_index] ) / ( delta_r * delta_r );
                break;
            }
            case BoundaryType::FixedDerivative :
            {
                return  ( ( previous_solution[boundary_index - 1] - previous_solution[boundary_index] )/delta_r  -  -1 * _fixed_temperature_derivative ) / ( delta_r );
                break;
            }
            default :
            {
                throw Error::UnknownBoundaryType;
            }
        }
    }
    //Fourth order spatial derivatives
    else if(order == ExplicitSolverSettings::SolverOrder::FOURTH)
    {
        //one away from right hand boundary
        if(index == -1)
        {
            
        }
        //one away from right hand boundary
        if(index == -2)
        {
            
        }
    }
}