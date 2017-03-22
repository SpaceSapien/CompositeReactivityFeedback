/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mixture.cpp
 * Author: chris
 * 
 * Created on March 20, 2017, 6:08 AM
 */

#include "Mixture.h"

using namespace MaterialLibrary;

Mixture::Mixture() 
{
    
}

bool Mixture::addCompound(const Compound &compound, const double &volume_fraction)
{
    _mixture.addComponent(compound,volume_fraction);
}

Composition<Compound> Mixture::getCompounds()
{
    return _mixture;
}

Composition<Isotope> Mixture::getIsotopeData()
{
    std::map<std::string,std::pair<Isotope,Real>> all_isotpes_composition;
    
    for( std::size_t compound_index = 0; compound_index < _mixture.size(); ++compound_index)
    {
        CompositionComponent<Compound> data = _mixture[compound_index];
        Compound current_compound = data._object;
        double volume_fraction = data._fraction;
        // Fixed Density at a temperature
        double density = current_compound.getDensity(500);
        double current_mass = volume_fraction * density;    
        double number_mols = current_mass / current_compound.getAtomicMass();
        
        std::vector<std::pair<Isotope,Real>> compound_mols = current_compound.getElementMols(number_mols);
        
        //Use a map to autocombine any duplicate elements
        for(std::size_t isotope_index = 0; isotope_index < compound_mols.size(); ++isotope_index)
        {
            Isotope current_isotope = compound_mols[isotope_index].first;
            
            std::string key = current_isotope.getZAID();
            
            if(all_isotpes_composition.find(key)  == all_isotpes_composition.end())
            {
                std::pair<Isotope,Real> current_isotope_data = compound_mols[isotope_index];
                all_isotpes_composition[key] = current_isotope_data;
            }
            else
            {
                all_isotpes_composition[key].second += compound_mols[isotope_index].second;
            }
        }
    }
    
    Composition<Isotope> homogenized_composition;
    
    for( auto isotope_data : all_isotpes_composition )
    {
        homogenized_composition.addComponent( isotope_data.second.first, isotope_data.second.second );
    }
    
    
    return homogenized_composition;
}

Mixture::~Mixture() 
{
    
}

