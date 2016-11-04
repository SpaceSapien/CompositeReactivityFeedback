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
    
    //Make sure the number of material equals the number of dimensions
    if( ! materials.size() == dimensions.size() )
    {
        throw Errors::IncorrectMaterialsDimensions;
    }
    
    this->_geometry = std::vector< std::pair<Materials,Dimension> >();
    
    
    // Check to make sure that the dimensions are increasing, if so create the 
    // material and dimensions vector
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
      
    return  MaterialLibrary::getMaterialProperties(material,T);
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