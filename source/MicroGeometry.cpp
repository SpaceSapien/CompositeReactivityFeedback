/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MicroGeometry.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 3:13 PM
 */
#include <iostream>
#include <cmath>
#include "MicroGeometry.h"

MicroGeometry::MicroGeometry(const std::vector<Materials> &materials,const std::vector<Dimension> &dimensions) 
{
    
    if( ! materials.size() == dimensions.size() )
    {
        throw Errors::IncorrectMaterialsDimensions;
    }
    
    this->_geometry = std::vector< std::pair<Materials,Dimension> >();
    
    Dimension last_dimension = 0;
    
    for( int index = 0; index < materials.size(); ++index )
    {
        if( ! dimensions[ index ] > last_dimension )
        {
            throw Errors::NonIncreasingDimensions;
        }
        last_dimension = dimensions[ index ];
        
        std::pair<Materials,Dimension> boundary = { materials[index], dimensions[ index] };
        this->_geometry.push_back( boundary );
    }
    
    this->_material_library = MaterialLibrary();    
}

MicroGeometry::MicroGeometry()
{
    
}

Dimension MicroGeometry::getOuterRadius() const
{
    return this->_geometry.back().second;
}

Dimension MicroGeometry::getFuelKernelRadius() const
{
    return this->_geometry.front().second;
}

void MicroGeometry::printGeometry()
{
    for( int index = 0; index < this->_geometry.size(); index++ )
    {
        std::cout<< index << "  " 
                 << getMaterialName( _geometry[index].first ) 
                 << "  " << _geometry[index].second << "\n";
    }
}

MaterialDataPacket MicroGeometry::getMaterialProperties(const Real &r,const Real &T)
{
    Materials material = this->getMaterial(r);
    
/*    if( material == MicroGeometry::_last_material && (T - 2) < MicroGeometry::_last_temperature && (T + 2) > MicroGeometry::_last_temperature )
    {
        return MicroGeometry::_last_packet;
    }
  */  
    
    std::pair<Real,Real> specific_heat_pair =this->_material_library.getSpecificHeatPair(material,T,0);
    //Fixed density ...
    std::pair<Real,Real> density_pair = this->_material_library.getDensityPair(material,1000,0); 
    std::pair<Real,Real> thermal_conductivity_pair = this->_material_library.getThermalConductivityPair(material,T,0);
    
    Real density = density_pair.first;
    Real specific_heat = specific_heat_pair.first;
    Real thermal_conductivity = thermal_conductivity_pair.first;
            
    Real thermal_conductivity_temperature_derivative = thermal_conductivity_pair.second;// this->_material_library.getThermalConductivityTemperatureDerivative(material,T,0);
    Real specific_heat_temperature_derivative = thermal_conductivity_pair.second; //this->_material_library.getSpecificHeatTemperatureDerivative(material,T,0);
    Real density_temperature_derivative = density_pair.second; //this->_material_library.getDensityTemperatureDerivative(material,T,0);
    
    MaterialDataPacket packet = MaterialDataPacket(thermal_conductivity, density, specific_heat, thermal_conductivity_temperature_derivative, specific_heat_temperature_derivative, density_temperature_derivative);
    
    //MicroGeometry::_last_temperature = T;
    //MicroGeometry::_last_material = material;
    //MicroGeometry::_last_packet = packet;
    
    return packet;
}



Materials MicroGeometry::getMaterial(const Real &r)
{
    Dimension location;
    int number_zones = _geometry.size();
    
    for(int index = 0; index < number_zones; index ++ )
    {
        if( r <= _geometry[index].second )
        {
            return _geometry[index].first;
        }
    }    
    
    throw Errors::rGreaterThanR;
    
}

Real MicroGeometry::getVolume()
{
    return 4.0/3.0 * M_PI * pow(this->getOuterRadius(),3);
}