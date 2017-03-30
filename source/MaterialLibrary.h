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
#include <map>
#include "EnumsAndFunctions.h"
#include "MaterialDataPacket.h"
#include "Composition.h"
#include "Isotope.h"
#include "Element.h"





namespace MaterialLibrary 
{
    
    typedef double Real;
    typedef double Dimension;
    

    class MicroGeometry;
    class Isotope;
    
    enum class Materials
    {
        U,
        UO2,
        DUO2,
        UN,
        UC,
        U3Si,
        SiC,
        C,
        Be,
        BeO,
        ZrB2,
        W,
        B4C,
        Mo,
        Nb,
        Zr,
        ZrO2,
        Graphene,
        Homogenous

    };
    
    enum class Atom
    {
        Be,
        B,
        C,
        Si,
        W,
        Nb,
        Mo,
        Zr,
        O,
        U,
        N
    };

    
    enum Errors
    {
        MaterialNotDefined,    
    };
    
    std::pair<Real,Real> getThermalConductivityPair(const Materials &material, const Real &T, const Real &dpa);
    std::pair<Real,Real> getDensityPair(const Materials &material,const Real &T,const Real &dpa);
    std::pair<Real,Real> getSpecificHeatPair(const Materials &material,const Real &T,const Real &dpa);
    
    void getMcnpMaterialCard(const Materials &material, const unsigned int &zone, const Real &average_temperature, std::string &material_card_entry, std::string &doppler_card,const Real &enrichment_fraction);
    std::map<std::string,bool>  getMcnpDopplerBroadenedCS(const std::vector<std::string> &cross_sections);
    void getMcnpMtMaterialCard(const Real &temperature,const std::vector<std::pair<int,Real>> &library_list,const std::string &library_base_name,std::string &library_name, Real &library_temperature);
    void getMcnpMTCard(const Materials &material, const Real &avg_temperature, std::string &mt_card_entry);
    std::map<std::string, std::pair<Real,std::string>> getMCNPCrossSectionLibraries(const Composition<Isotope> &pre_modified_libraries);
    
    Real getMeltingPoint(const Materials &material);
    Real getLinearExpansionCoeifficient(const Materials &material,const Real &T,const Real &dpa);
    MaterialDataPacket getMaterialProperties(const Materials &material,const Real &T);
    Real interpolateDataArrays(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);
    std::pair<Real,Real> interpolateDataAndTemperatureArraysAndDerivative(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x);
    
    std::string getMaterialName(const Materials &material);
    Materials getMaterialFromName(const std::string &material_name);
    
    

};



#endif /* MATERIALIBRARY_H */

