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

class InfiniteCompositeReactor;

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
    
    InfiniteCompositeReactor* _reactor;
    std::vector<Real> _delayed_precursors;
    Real _initial_power;
    Real _current_power;
    DelayedNeutronSet _delayed_neutron_set;
    const Real _power_per_fission = 3.20435e-11; //Joules per fission
    Real _current_time = 0;
    Real _kinetics_time_step;
    
    ReactorKinetics();
    ReactorKinetics(InfiniteCompositeReactor* reactor,const Real &initial_power, const DelayedPrecursorInitialState &state );
    Real solveForPower( const Real &simulation_time);
    
private:
    DelayedNeutronSet getDelayedNeutronInfo(const std::vector< std::pair<FissionableIsotope,Real> > &fission_listing);
    

};

#endif /* REACTORKINETICS_H */

