/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaterialDataPacket.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 5:08 PM
 */
#include <iostream>
#include "MaterialDataPacket.h"

MaterialDataPacket::MaterialDataPacket(const Real thermal_conductivity,const Real density, const Real specific_heat, 
        const Real thermal_conductivity_temperature_derivative,const Real specific_heat_temperature_derivative, const Real density_temperature_derivative)
{
    _thermal_conductivity = thermal_conductivity;
    _density = density;
    _specific_heat = specific_heat;
    _thermal_conductivity_temperature_derivative = thermal_conductivity_temperature_derivative;
    _specific_heat_temperature_derivative = specific_heat_temperature_derivative;
    _density_temperature_derivative = density_temperature_derivative;
}

void MaterialDataPacket::printDataPacket()
{
    std::cout<< "Thermal Conductivity : " << this->_thermal_conductivity << " W/m-K\n";
    std::cout<< "Density : " << this->_density << " kg/m^3\n";
    std::cout<< "Specific Heat : " << this->_specific_heat << " J/kg-K\n";
    std::cout<< "Specific Heat Temp Derivative : " << this->_specific_heat_temperature_derivative << " J/m-(K^2)\n";
    std::cout<< "Thermal Conductivity Temp Derivative : " << this->_thermal_conductivity_temperature_derivative << " W/m-(K^2)\n";
    std::cout<< "Density Temp Derivative : " << this->_density_temperature_derivative << " kg/m^3-K\n";
}



