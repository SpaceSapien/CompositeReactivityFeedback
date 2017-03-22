/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MicroGeometry.h
 * Author: chris
 *
 * Created on November 20, 2015, 3:13 PM
 */



#ifndef MICROGEOMETRY_H
#define MICROGEOMETRY_H
#include <vector>
#include "EnumsAndFunctions.h"
#include "MaterialLibrary.h"
#include "MaterialDataPacket.h"

namespace MaterialLibrary
{

class MicroGeometry 
{

public:
    
    enum Errors
    {
        IncorrectMaterialsDimensions,
        NonIncreasingDimensions,
        rGreaterThanR
    };
    
    
    MicroGeometry(const std::vector<Materials> &materials, const std::vector<Dimension> &geometry);
    MicroGeometry();
    
    std::pair<Real,Real> getSpecificHeatPair(const Real &r,const  Real &T);
    std::pair<Real,Real> getDensityPair(const Real &r,const  Real &T);
    std::pair<Real,Real> getThermalConductivityPair(const Real &r,const  Real &T);
    
    
    MaterialDataPacket getMaterialProperties(const Real &r,const Real &T);
    Materials getMaterial(const Real &r);
    Dimension getOuterRadius() const;
    Dimension getFuelKernelRadius() const;
    std::vector<Real> getVolumeFractions() const;
    void printGeometry();
    
    MaterialDataPacket getHomogenizedMaterialProperties(const Real &T);
    void getHomogenizedMcnpMaterialCard(const int &zone, const std::vector<Real> &geometry_temperature_data, const Real &enrichment_fraction,
           std::string &material_cards,std::string &doppler_cards, std::string &mt_cards) const;
    
    
    static MaterialDataPacket _last_packet;
    static Real _last_temperature;
    static Materials _last_material;
    
    std::vector< std::pair<Materials,Dimension> > _geometry;
    
private:
    
    
    

};

}

#endif /* MICROGEOMETRY_H */

