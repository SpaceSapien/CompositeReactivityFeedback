/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinMonteCarlo.h
 * Author: chris
 *
 * Created on April 5, 2017, 3:24 AM
 */

#ifndef FUELPINMONTECARLO_H
#define FUELPINMONTECARLO_H

#include <string>
#include "ReactorMonteCarlo.h"
#include "FuelPinReactor.h"

class FuelPinReactor;

class FuelPinMonteCarlo : public ReactorMonteCarlo
{
public:
    
    FuelPinMonteCarlo(FuelPinReactor* reactor, const std::string &run_dir);    
    virtual ~FuelPinMonteCarlo();
    
    Real _coolant_channel_radius;
    Real _fuel_pin_pitch;
    Real _height;
    int _number_macro_cells;
    
    virtual std::string getMaterialCards();
    virtual std::string getCellCards();
    virtual std::string getSurfaceCards();    
    //virtual std::vector< std::vector<Real> > getZoneCellRelativePowerDensity();
    
private:
    
    FuelPinReactor* _reactor;

};

#endif /* FUELPINMONTECARLO_H */

