/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MicroCellBoundaryCondition.h
 * Author: chris
 *
 * Created on December 17, 2015, 4:22 PM
 */

#ifndef MICROCELLBOUNDARYCONDITION_H
#define MICROCELLBOUNDARYCONDITION_H

#include "EnumsAndFunctions.h"
#include "MicroCell.h"
class MicroCell;



class MicroCellBoundaryCondition 
{
    
public:
    
    enum BoundaryType
    {
        FixedTemperature,
        ReflectedHeatFlux,
        FixedHeatFlux
    };
    
    enum Error
    {
        UnknownBoundaryType        
    };

    
    
    BoundaryType _boundary;
    
    static MicroCellBoundaryCondition* getReflectedHeatFluxBoundaryConditionFactory();
    static MicroCellBoundaryCondition* getFixedTemperatureBoundaryConditionFactory(const Real &fixed_temperature);
    static MicroCellBoundaryCondition* getFixedHeatFluxBoundaryConditionFactory(const Real &fixed_heat_flux);
    
    MicroCellBoundaryCondition();
    
    Real getHeatFlux(MicroCell * microcell, const Real &inward_heat_flux = 0) const;

private:
    
    MicroCellBoundaryCondition(const BoundaryType &type);
    Real _fixed_temperature;
    Real _fixed_heat_flux;
};



#endif /* MICROCELLBOUNDARYCONDITION_H */

