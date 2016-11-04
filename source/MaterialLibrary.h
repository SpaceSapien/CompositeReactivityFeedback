/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Materialibrary.h
 * Author: chris
 *
 * Created on November 18, 2015, 2:02 PM
 */

#ifndef MATERIALIBRARY_H
#define MATERIALIBRARY_H
#include <vector>
#include "EnumsAndFunctions.h"
#include "MaterialDataPacket.h"

class MaterialLibrary 
{
public:
    
    enum Errors
    {
        MaterialNotDefined,    
    };
    
    MaterialLibrary();
    
    std::pair<Real,Real> static getThermalConductivityPair(const Materials &material, const Real &T, const Real &dpa);
    std::pair<Real,Real> static getDensityPair(const Materials &material,const Real &T,const Real &dpa);
    std::pair<Real,Real> static getSpecificHeatPair(const Materials &material,const Real &T,const Real &dpa);
    
    void static getMcnpMaterialCard(const Materials &material, const unsigned int &zone, std::string &material_card_entry, std::string &doppler_card,const Real &enrichment_fraction);
 
    
    Real static getMeltingPoint(const Materials &material);
    Real static getLinearExpansionCoeifficient(const Materials &material,const Real &T,const Real &dpa);
    MaterialDataPacket static getMaterialProperties(const Materials &material,const Real &T);
    
   
    
private:
    
    Real static interpolateDataArrays(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);
    std::pair<Real,Real> static interpolateDataAndTemperatureArraysAndDerivative(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);
    

};

#endif /* MATERIALIBRARY_H */

