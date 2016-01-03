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
#include "ExplicitSolverSettings.h"


class MicroCellBoundaryCondition 
{
    
public:
    
    enum BoundaryType
    {
        Fixed,
        Reflected,
        FixedDerivative
    };
    
    enum Error
    {
        UnknownBoundaryType        
    };

    
    BoundaryType _boundary;
    
    MicroCellBoundaryCondition static getRefelectedBoundaryCondition();
    MicroCellBoundaryCondition static getFixedBoundaryCondition(const Real &fixed_temperature);
    MicroCellBoundaryCondition static getFixedDerivativeBoundaryCondition(const Real &fixed_spatial_derivative);
    MicroCellBoundaryCondition();
    
    Real getdTdr(const ExplicitSolverSettings::SolverOrder &order,const int &index,const std::vector<Real> &previous_solution, const std::vector<Real> &grid );
    Real getd2Tdr2(const ExplicitSolverSettings::SolverOrder &order,const int &index,const std::vector<Real> &previous_solution, const std::vector<Real> &grid );

private:
    
    MicroCellBoundaryCondition(const BoundaryType &type);
    Real _fixed_temperature;
    Real _fixed_temperature_derivative;
};



#endif /* MICROCELLBOUNDARYCONDITION_H */

