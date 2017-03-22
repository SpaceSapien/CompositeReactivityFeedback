/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Compound.h
 * Author: chris
 *
 * Created on March 20, 2017, 4:33 AM
 */

#ifndef COMPOUND_H
#define COMPOUND_H
#include "Composition.h"
#include "Element.h"
#include "MaterialLibrary.h"

namespace MaterialLibrary
{

    class Compound
    {
    public:

        Compound(const Materials &material);
        ~Compound();

        bool addElement(const Element &element,const Real &relative_atomic_density);      
        Real getAtomicMass();
        static Compound getStandardCompound(const Materials &material, const double &enrichment);
        std::vector<std::pair<Isotope,Real>> getElementMols(const Real &compound_mols);       
        
        Materials getMaterial();
        Real getDensity(const Real &T) const;
        
    private:

        Composition<Element> _element_mol_fractions;
        std::string _compound;
        Materials _material;
        
    };
}

#endif /* COMPOUND_H */

