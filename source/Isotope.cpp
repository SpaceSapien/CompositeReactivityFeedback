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
    
    //Carbon 
    atomic_mass["6012"] = 12.0;
    atomic_mass["6013"] = 13.003355;
    //Oxygen
    atomic_mass["8016"] = 15.994915;
    atomic_mass["8017"] = 16.999132;
    atomic_mass["8018"] = 17.999160;
    // Zirconium
    atomic_mass["40090"] = 89.904704;
    atomic_mass["40091"] = 90.905645;
    atomic_mass["40092"] = 91.905040;
    atomic_mass["40094"] = 93.906316;
    atomic_mass["40096"] = 95.908276;
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