/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinCompositeMicroCell.h
 * Author: chris
 *
 * Created on April 6, 2017, 6:18 PM
 */

#ifndef FUELPINCOMPOSITEMICROCELL_H
#define FUELPINCOMPOSITEMICROCELL_H
#include "CompositeMicroCell.h"
#include "FuelPinCompositeMicroCell.h"
#include "FuelPinReactor.h"
#include "InputDataFunctions.h"

class FuelPinReactor;

class FuelPinCompositeMicroCell : public CompositeMicroCell
{
public:
    
    Real _surface_resistance;
    Real _coating_surface_area;
    Real _number_particles_per_volume;
    Real _matrix_volume_fraction;
    Real _macro_scale_position;
    
    FuelPinCompositeMicroCell(FuelPinReactor* reactor, const Real &temperature,const Real &macroscale_position);
    virtual ~FuelPinCompositeMicroCell();
    virtual std::vector<MicroSolution> iterateInitialConditions(const std::vector<Real> &power_distribution);
    
    Real solveMicroCellCoupling(const Real &current_macrocell_temperature, const Real &microcell_internal_power_density, const Real &time_step);
    
private:
    
    virtual void solveAndSetSurfaceResistance(const Real &cell_power);

};

#endif /* FUELPINCOMPOSITEMICROCELL_H */

