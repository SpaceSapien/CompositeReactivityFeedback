/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Element.cpp
 * Author: chris
 * 
 * Created on March 20, 2017, 2:33 AM
 */

#include "Element.h"

using namespace MaterialLibrary;

Element::Element(const std::string &name,const std::string &symbol,const unsigned int atomic_number, const Atom &atom)
{
    _element_name = name;
    _atomic_number = atomic_number;
    _symbol = symbol;
    _atom = atom;
    
}

Element::~Element() {}

double Element::getAtomicMass() const
{
    double average_atomic_mass;

    for( std::size_t index = 0; index < _isotope_composition.size(); ++index)
    {
        
        CompositionComponent<Isotope> isotope_composition = _isotope_composition[index];
        average_atomic_mass += isotope_composition._fraction * isotope_composition._object.getAtomicMass();
    }

    return average_atomic_mass;
}

std::string Element::getSymbol() const
{
    return _symbol;
}

Atom Element::getAtom() const
{
    return _atom;
}

int Element::getAtomicNumber() const
{
    return _atomic_number;
}



std::string Element::getElementName() const
{
    return _element_name;
}

Element Element::getEnrichedUranium(const double &enrichment)
{
    std::string name = "Uranium";
    std::string symbol = "U";
    int atomic_number = 92;
    
    if(enrichment > 1 || enrichment < 0)
    {
        throw std::string("Enrichment fraction must be between 0 and 1");
    }
    
    std::vector<int> nucleon_numbers = { 234,      235,      238      };
    std::vector<double> fractions    = { 0.00008 * enrichment, enrichment,   1 - 1.008 * enrichment };
    
    //Create the natural elements based on the data
    Element element = Element(name,symbol,atomic_number, Atom::U); 
    
    for(std::size_t index=0; index < nucleon_numbers.size(); ++index)
    {
        Isotope isotope = Isotope(&element, nucleon_numbers[index]);
        element.addIsotope(isotope, fractions[index]);
    }
    
    return element;
}

std::vector<std::pair<Isotope,Real>> Element::getIsotopeMols(const Real &element_mols ) const
{
    std::vector<std::pair<Isotope,Real>> isotope_data; 
    Real density = element_mols * getAtomicMass();
    
    for(std::size_t isotope_index=0; isotope_index < _isotope_composition.size(); ++isotope_index)
    {
        CompositionComponent<Isotope> current_isotope = _isotope_composition[isotope_index];
        Real isotope_density = current_isotope._fraction * density;
        Real isotope_mols = isotope_density / current_isotope._object.getAtomicMass();
        std::pair<Isotope,Real> current_data = {current_isotope._object, isotope_mols};
        isotope_data.push_back(current_data);
    }
    
    return isotope_data;
}

Element Element::getNaturalAtomicData(Atom atom)
{
    std::vector<int> nucleon_numbers;
    std::vector<double> fractions;
    
    std::vector<int> mcnp_nucleons;
    std::vector<double> mcnp_fractions;
    
    
    std::string name;
    std::string symbol;
    int atomic_number;
    
  
    
    switch(atom)
    {
        case Atom::C :
        {
            name = "Carbon";
            symbol = "C";
            atomic_number = 6;

            nucleon_numbers = { 12,     13     };
            fractions =       { 0.9893, 0.0107 };          
            
            
            break;
        }
        case Atom::O :
        {
            name = "Oxygen";
            symbol = "O";
            atomic_number = 8;
           
            nucleon_numbers = { 16,      17,      18      };
            fractions =       { 0.99757, 0.00038, 0.00205 };
            
            break;
        }
        
        case Atom::U :
        {
            name = "Uranium";
            symbol = "U";
            atomic_number = 92;

            nucleon_numbers = { 234,      235,      238      };
            fractions =       { 0.000055, 0.0072,   0.992745 };
            
            break;
        }
        
        case Atom::Zr :
        {
            name = "Zirconium";
            symbol = "Zr";
            atomic_number = 40;

            nucleon_numbers = { 90,      91,      92,      94 ,    96 };
            fractions =       { 0.5145,  0.1122,  0.1715,  0.1738, 0.0280};
            
            break;
        }
    }
    
    //Create the natural elements based on the data
    Element element = Element(name,symbol,atomic_number, atom); 
    
    for(std::size_t index=0; index < nucleon_numbers.size(); ++index)
    {
        Isotope isotope = Isotope(&element, nucleon_numbers[index]);
        element.addIsotope(isotope, fractions[index]);
    }
    
    return element;
}


bool Element::addIsotope(const Isotope &isotope,const double &fraction)
{
    _isotope_composition.addComponent(isotope,fraction);
}

    
