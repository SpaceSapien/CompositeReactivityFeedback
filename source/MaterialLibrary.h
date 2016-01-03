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

class MaterialLibrary 
{
public:
    
    enum Errors
    {
        MaterialNotDefined,    
    };
    
    MaterialLibrary();
    MaterialLibrary(const MaterialLibrary& orig);
    
    std::pair<Real,Real> getThermalConductivityPair(const Materials &material, const Real &T, const Real &dpa);
    std::pair<Real,Real> getDensityPair(const Materials &material,const Real &T,const Real &dpa);
    std::pair<Real,Real> getSpecificHeatPair(const Materials &material,const Real &T,const Real &dpa);
    
    Real getMeltingPoint(const Materials &material);
    Real getLinearExpansionCoeifficient(const Materials &material,const Real &T,const Real &dpa);
    
   
    
private:
    
    Real interpolateDataArrays(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);
    std::pair<Real,Real> interpolateDataAndTemperatureArraysAndDerivative(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);


};

#endif /* MATERIALIBRARY_H */

