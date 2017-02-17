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


MicroCellBoundaryCondition*  MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory()
{
    MicroCellBoundaryCondition* boundary_condition = new MicroCellBoundaryCondition(BoundaryType::ReflectedHeatFlux);
    return boundary_condition; 
}

/**
 * Factory function that makes a new boundary conidtion
 * 
 * @param fixed_spatial_derivative dr/dT = cost where const is fixed_spatial_derivative
 * @return MicroCellBoundaryCondition
 */
MicroCellBoundaryCondition*  MicroCellBoundaryCondition::getFixedHeatFluxBoundaryConditionFactory(const Real &fixed_heat_flux)
{
    MicroCellBoundaryCondition* boundary_condition = new MicroCellBoundaryCondition(BoundaryType::FixedHeatFlux);
    boundary_condition->_fixed_heat_flux = fixed_heat_flux;
    return boundary_condition; 
}

    
MicroCellBoundaryCondition* MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(const Real &fixed_temperature)
{
    MicroCellBoundaryCondition* boundary_condition = new MicroCellBoundaryCondition(BoundaryType::FixedTemperature);
    boundary_condition->_fixed_temperature = fixed_temperature;
    return boundary_condition; 
}

Real MicroCellBoundaryCondition::getHeatFlux(MicroCell*  micro_cell, const Real &inward_heat_flux ) const
{
    switch(_boundary)
    {
        case BoundaryType::FixedHeatFlux :
        {
            return  -1 * this->_fixed_heat_flux;
            break;
        }
        case BoundaryType::ReflectedHeatFlux :
        {
            return  0;
            break;
        }
        case BoundaryType::FixedTemperature :
        {
            return  -inward_heat_flux;
            break;
        }
        default :
        {
            throw Error::UnknownBoundaryType;
        }
    } 
}

/**
 * 
 * @param boundary_type FixedTemperature, ReflectexdHeatFlux, or FixedHeatFlux
 * @return the boundary type
 */
MicroCellBoundaryCondition::BoundaryType MicroCellBoundaryCondition::getBoudaryTypeFromString(const std::string &boundary_type)
{
    if(boundary_type == "FixedTemperature")
    {
        return MicroCellBoundaryCondition::FixedTemperature;
    }
    if(boundary_type == "ReflectedHeatFlux")
    {
        return MicroCellBoundaryCondition::ReflectedHeatFlux;
    }
    if(boundary_type == "FixedHeatFlux")
    {
        return MicroCellBoundaryCondition::FixedHeatFlux;
    }
    
    throw "Unknown Boundary Condition: " + boundary_type + " line " + std::to_string(__LINE__) + " file " + __FILE__;
}