/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Isotope.cpp
 * Author: chris
 * 
 * Created on March 20, 2017, 2:24 AM
 */

#include "Isotope.h"
#include "Element.h"
#include <map>

using namespace MaterialLibrary;

Isotope::Isotope(){}

Isotope::Isotope(Element* element,const int nucleon_number)
{
    _nucleon_number = nucleon_number;
    _atomic_number = element->getAtomicNumber();
    _atom = element->getAtom();
    _element_name = element->getElementName();
    _element_symbol = element->getSymbol();
    
    std::string zaid = this->getZAID();
    
    _atomic_mass = Isotope::getAtomicMassFomZAID(zaid);
}

std::string Isotope::getElementSymbol() const
{
    return _element_symbol;
}


std::string Isotope::getElementName() const
{
    return _element_name;
}

Atom Isotope::getAtom() const
{
    return _atom;
}

int Isotope::getAtomicNumber() const
{
    return _atomic_number;
}

int Isotope::getNucleonNumber() const
{
    return _nucleon_number;
}

double Isotope::getAtomicMass() const
{
    return _atomic_mass;
}

std::string Isotope::getZAID() const
{
    std::string atomic_mass_str = std::to_string( this->getNucleonNumber() );
    int string_length = atomic_mass_str.length();

    for(int current_length = string_length; current_length < 3; ++current_length)
    {
        atomic_mass_str = "0" + atomic_mass_str;
    }

    int atomic_number = this->getAtomicNumber();
    
    std::string zaid_string = std::to_string( atomic_number) + atomic_mass_str;
    return zaid_string;
}

double Isotope::getAtomicMassFomZAID(const std::string &zaid)
{
    std::map<std::string,double> atomic_mass;
    
    //Boron
    atomic_mass["5010"] = 10.012937;
    atomic_mass["5011"] = 11.009305;
    
    //Beryllium
    atomic_mass["9004"] = 9.012182;
    
    //Carbon 
    atomic_mass["6012"] = 12.0;
    atomic_mass["6013"] = 13.003355;
    
    //Nitrogen
    atomic_mass["7014"] = 14.003074;
    atomic_mass["7015"] = 15.000109;
    
    //Oxygen
    atomic_mass["8016"] = 15.994915;
    atomic_mass["8017"] = 16.999132;
    atomic_mass["8018"] = 17.999160;
    
    //Silicon
    atomic_mass["14028"] = 27.976927;
    atomic_mass["14029"] = 28.976495;
    atomic_mass["14030"] = 29.973770;
    
    // Zirconium
    atomic_mass["40090"] = 89.904704;
    atomic_mass["40091"] = 90.905645;
    atomic_mass["40092"] = 91.905040;
    atomic_mass["40094"] = 93.906316;
    atomic_mass["40096"] = 95.908276;
    
    // Niobium
    atomic_mass["41093"] = 92.906378;
    
    // Molybdenum
    atomic_mass["42092"] = 91.906810;
    atomic_mass["42094"] = 93.905088;
    atomic_mass["42095"] = 94.905841;
    atomic_mass["42096"] = 95.904679;
    atomic_mass["42097"] = 96.906021;
    atomic_mass["42098"] = 97.905408;
    atomic_mass["42100"] = 99.907477;
    
    // Tungsten
    atomic_mass["74180"] = 179.946706;
    atomic_mass["74182"] = 181.948206;
    atomic_mass["74183"] = 182.950224;
    atomic_mass["74184"] = 183.950933;
    atomic_mass["74186"] = 185.954362;
    
    //Uranium
    atomic_mass["92234"] = 234.040946;
    atomic_mass["92235"] = 235.043923;
    atomic_mass["92238"] = 238.050783;
    
    if ( atomic_mass.find(zaid) == atomic_mass.end() ) 
    {
        throw "ZAID " + zaid + " not defined " + std::to_string(__LINE__) + " " + __FILE__;
    } 
    else 
    {
        return atomic_mass[zaid];
    }
    
    
}