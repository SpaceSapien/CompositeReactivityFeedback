/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Element.h
 * Author: chris
 *
 * Created on March 20, 2017, 2:33 AM
 */

#ifndef ELEMENT_H
#define ELEMENT_H
#include "Isotope.h"
#include "Composition.h"
#include "MaterialLibrary.h"



namespace MaterialLibrary
{
    enum class Atom;
    class Isotope;
    
    class Element
    {
        

    public:
        
        

        //Constructors
        Element(const std::string &name,const std::string &symbol,const unsigned int atomic_number,const Atom &atom);
        virtual ~Element();
        //Getters
        double getAtomicMass() const;
        int getAtomicNumber() const;
        std::string getElementName() const;   
        std::string getSymbol() const;
        Atom getAtom() const;
        
        std::vector<std::pair<Isotope,double>> getIsotopeMols(const double &element_mols) const;
        //Add
        bool addIsotope(const Isotope &isotope,const double &fraction);
        //Factories
        static Element getNaturalAtomicData(Atom atom);
        static Element getEnrichedUranium(const double &enrichment);

    private:
        Composition<Isotope> _isotope_composition;  
        std::string _element_name;
        std::string _symbol;
        int _atomic_number;
        Atom _atom;

    };
    
    
    
}


#endif /* ELEMENT_H */

