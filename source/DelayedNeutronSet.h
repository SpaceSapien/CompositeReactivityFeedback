/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DelayedNeutronSet.h
 * Author: chris
 *
 * Created on December 8, 2015, 8:04 PM
 */


#ifndef DELAYEDNEUTRONSET_H
#define DELAYEDNEUTRONSET_H
#include "EnumsAndFunctions.h"
#include <vector>
#include <string>

class DelayedNeutronSet 
{
    
public:
    
    enum Error
    {
        UnknownIsotope,
        NonMatchingDelayedTimeConstantsAndFractions
    };
    
    const std::vector<Real> _beta_time_constants= { 0.0129, 0.0311, 0.134, 0.331, 1.26, 3.21 };
    
    std::vector<Real> _delayed_neutrons_per_fission_groups;
    Real _delayed_neutrons_per_fission;
    Real _fast_neutrons_per_fission;
    Real _neutrons_per_fission;
    
    FissionableIsotope _isotope_name;
    
    DelayedNeutronSet();
    DelayedNeutronSet(const std::vector< std::pair<FissionableIsotope,Real>> &fission_talley);
    DelayedNeutronSet(const FissionableIsotope &isotope_name);
    DelayedNeutronSet operator=(const DelayedNeutronSet &isotope_name);   
    Real getTotalBeta();
    int numberDelayedNeutronSets();
    void print();
    
private:
    
    Real sumVector(std::vector<Real> data);

};

#endif /* DELAYEDNEUTRONSET_H */

