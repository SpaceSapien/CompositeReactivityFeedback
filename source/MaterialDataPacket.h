/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaterialDataPacket.h
 * Author: chris
 *
 * Created on November 20, 2015, 5:08 PM
 */


#ifndef MATERIALDATAPACKET_H
#define MATERIALDATAPACKET_H
#include "EnumsAndFunctions.h"
class MaterialDataPacket {
public:
    MaterialDataPacket(const Real thermal_conductivity,const Real density, const Real specific_heat, 
        const Real thermal_conductivity_temperature_derivative,const Real specific_heat_temperature_derivative, const Real density_temperature_derivative);
        
    void printDataPacket();
    
    
    Real _thermal_conductivity;
    Real _density;
    Real _specific_heat;
    Real _thermal_conductivity_temperature_derivative;
    Real _specific_heat_temperature_derivative;
    Real _density_temperature_derivative;
       
    
    
    
private:

};

#endif /* MATERIALDATAPACKET_H */

