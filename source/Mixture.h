/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mixture.h
 * Author: chris
 *
 * Created on March 20, 2017, 6:08 AM
 */

#ifndef MIXTURE_H
#define MIXTURE_H
#include "MaterialLibrary.h"
#include "Compound.h"
#include "Isotope.h"

namespace MaterialLibrary
{

    class Mixture
    {
        public:
            Mixture();
            bool addCompound(const Compound &compound, const double &volume_fraction);
            ~Mixture(); 
            Composition<Isotope> getIsotopeData();
            Composition<Compound> getCompounds();
            
        private:
            
            Composition<Compound> _mixture;
    };

}
#endif /* MIXTURE_H */

