/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Isotope.h
 * Author: chris
 *
 * Created on March 20, 2017, 2:24 AM
 */

#ifndef ISOTOPE_H
#define ISOTOPE_H
#include <string>


#include "MaterialLibrary.h"

namespace MaterialLibrary
{
    class Element;    
    enum class Atom;
    
    class Isotope
    {

    public:

        Isotope();
        Isotope(Element* element,const int nucleon_number);    
        double getAtomicMass() const;
        std::string getZAID() const;
        Atom getAtom() const;
        std::string getElementName() const;
        int getAtomicNumber() const;
        int getNucleonNumber() const;
        std::string getElementSymbol() const;
        
        
        static double getAtomicMassFomZAID(const std::string &zaid); 

    private:

        std::string _element_name;
        std::string _element_symbol;
        
        double _atomic_mass;
        int _nucleon_number;
        int _atomic_number;
        Atom _atom; 

    };
}
#endif /* ISOTOPE_H */

