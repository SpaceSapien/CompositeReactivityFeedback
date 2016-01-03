/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactorKinetics.h
 * Author: chris
 *
 * Created on December 8, 2015, 3:50 PM
 */

#ifndef REACTORKINETICS_H
#define REACTORKINETICS_H
#include "EnumsAndFunctions.h"
#include "DelayedNeutronSet.h"

class ReactorKinetics 
{
public:
   
    enum Error
    {
        UnknownDelayedPrecoursorInitialState
    };
    
    enum DelayedPrecursorInitialState
    {
        NoInitialPrecursors,
        EquilibriumPrecursors
    };
    
    
    std::vector<Real> _delayed_precursors;
    Real _initial_power;
    Real _current_power;
    DelayedNeutronSet _delayed_neutron_set;
    const Real _power_per_fission = 3.20435e-11; //Joules per fission
    Real _current_time = 0;
    
    ReactorKinetics();
    ReactorKinetics(const Real &initial_power, const DelayedPrecursorInitialState &state );
    Real solveForPower
    (
        const Real &simulation_time, 
        const Real &k_effective, 
        const Real &neutron_generation_time, 
        const std::vector< std::pair<FissionableIsotope,Real> > &fission_listing,
        std::vector< std::pair<Real,Real>> &power_record,
        std::vector<std::pair<Real,std::vector<Real>>> &delayed_record
    );
    
private:
    DelayedNeutronSet getDelayedNeutronInfo(const std::vector< std::pair<FissionableIsotope,Real> > &fission_listing);
    

};

#endif /* REACTORKINETICS_H */

