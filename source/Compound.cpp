/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Compound.cpp
 * Author: chris
 * 
 * Created on March 20, 2017, 4:33 AM
 */
#include "Compound.h"
#include "Element.h"
#include "MaterialLibrary.h"

using namespace MaterialLibrary;

Compound::Compound(const Materials &material) 
{
    _material = material;
    std::string compound_name = MaterialLibrary::getMaterialName(material);    
    _compound = compound_name;
}

bool Compound::addElement(const Element &element,const double &relative_atomic_density)
{
    _element_mol_fractions.addComponent(element,relative_atomic_density);
}

double Compound::getDensity(const double &T) const
{
    std::pair<double,double> density = MaterialLibrary::getDensityPair(this->_material,T,0);
    return density.first;
}


Compound Compound::getStandardCompound(const Materials &material, const double &enrichment)
{
    Compound compound = Compound(material);    
    std::vector<Element> elemental_composition;
    std::vector<double> amounts;
    
    switch(material)
    {
        case Materials::C :
        {
            Element carbon = Element::getNaturalAtomicData(Atom::C);
            elemental_composition = { carbon };
            amounts = { 1.0 };
            break;
        }
        case Materials::ZrO2 :
        {
            Element oxygen = Element::getNaturalAtomicData(Atom::O);
            Element zirconium = Element::getNaturalAtomicData(Atom::Zr);            
            
            elemental_composition = { zirconium, oxygen };
            amounts = { 1.0, 2.0 };
            
            break;
        }
        case Materials::UO2 :
        {
            Element oxygen = Element::getNaturalAtomicData(Atom::O);
            Element uranium = Element::getEnrichedUranium(enrichment);            
            
            elemental_composition = { uranium, oxygen };
            amounts = { 1.0, 2.0 };
            
            break;
        }
    }
    
    for(std::size_t index=0; index < elemental_composition.size(); ++index)
    {
        Element current_element = elemental_composition[index];
        compound.addElement(current_element, amounts[index]);
    }
    
    return compound;
}

Materials Compound::getMaterial()
{
    return _material;
}

Compound::~Compound() {}


Real Compound::getAtomicMass()
{
    std::size_t number_elements = this->_element_mol_fractions.size();
    Real total_atomic_mass = 0;
    
    for(std::size_t index = 0; index < number_elements; ++index )
    {
        CompositionComponent<Element> data = _element_mol_fractions[index];
        total_atomic_mass += data._object.getAtomicMass() * data._amount;        
    }
    
    return total_atomic_mass;
}


std::vector<std::pair<Isotope,Real>> Compound::getElementMols(const Real &compound_mols)
{
    std::vector<std::pair<Isotope,Real>> isotope_mols;
            
    std::size_t number_elements = _element_mol_fractions.size();
    
    for( std::size_t index = 0; index < number_elements; ++index)
    {
        CompositionComponent<Element> data = _element_mol_fractions[index];
        Element current_element = data._object;
        Real element_mols = compound_mols * data._amount;
        
        std::vector<std::pair<Isotope,Real>> current_element_isotope_mols = current_element.getIsotopeMols(element_mols);
        
        for(std::size_t add_index = 0; add_index < current_element_isotope_mols.size(); add_index++)
        {
            isotope_mols.push_back(current_element_isotope_mols[add_index]);
        }
    }
    return isotope_mols;
}