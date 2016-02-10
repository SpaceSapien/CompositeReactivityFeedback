/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactorKinetics.cpp
 * Author: chris
 * 
 * Created on December 8, 2015, 3:50 PM
 */
#include "DelayedNeutronSet.h"
#include "ReactorKinetics.h"
#include <iostream>

ReactorKinetics::ReactorKinetics() { }


/**
 * Constructor to get the initial conditions for the reactor kinetics equations
 * @param initial_power
 * @param state
 */
ReactorKinetics::ReactorKinetics(const Real &initial_power,const DelayedPrecursorInitialState &state) 
{
    _initial_power = initial_power;
    _current_power = _initial_power;
    _current_time  = 0;
    
    //default delayed neutron set
    _delayed_neutron_set = DelayedNeutronSet(FissionableIsotope::U235);
    
    
    switch(state)
    {
        case DelayedPrecursorInitialState::NoInitialPrecursors :
        {
            //All delayed precursors are set to zero
            _delayed_precursors.resize(_delayed_neutron_set.numberDelayedNeutronSets(),0);
            break;
        }
        case DelayedPrecursorInitialState::EquilibriumPrecursors :
        {
            //Delayed precursors are set to their equilibrium states as if the reactor
            //operated for a long time at that power
            _delayed_precursors.resize(_delayed_neutron_set.numberDelayedNeutronSets(),0);
            
            //initialize the delated precoursors
            for(int index = 0; index < _delayed_precursors.size(); index++)
            {
                Real nu_dk = _delayed_neutron_set._delayed_neutrons_per_fission_groups[index];
                Real lambda_k = _delayed_neutron_set._beta_time_constants[index];
                _delayed_precursors[index] = (_initial_power * nu_dk) / ( _power_per_fission * lambda_k );
            }
            
            break;
        }
        default :
        {
            throw UnknownDelayedPrecoursorInitialState;
        }
    }
    
       
}

/**
 * Solve the reactor kinetics equations over a given time period
 * 
 * @param simulation_coupled_time_step 
 * @param simulation_time_step
 * @param k_effective
 * @param neutron_generation_time
 * @param fission_listing
 * @param power_record
 * @param time_record
 * @param delayed_record
 * @return 
 */
Real ReactorKinetics::solveForPower
(
    const Real &simulation_coupled_time_step, 
    const Real &k_effective, 
    const Real &neutron_generation_time, 
    const std::vector< std::pair<FissionableIsotope,Real> > &fission_listing,
    std::vector< std::pair<Real,Real>> &power_record,
    std::vector<std::pair<Real,std::vector<Real>>> &delayed_record
)
{
    //Get the proper beta based on the fissions
    _delayed_neutron_set = DelayedNeutronSet(fission_listing);
    //_delayed_neutron_set.print();
    
    //Calculate some parameters
    Real reactivity = ( k_effective - 1 )/k_effective;
    Real total_beta = _delayed_neutron_set.getTotalBeta();
    Real total_neutrons_per_fission = _delayed_neutron_set._neutrons_per_fission;
    
    long int iterations = 0;
    
    Real simulation_time_step = 10e-9;
    //Loop over a simulation coupled time step
    Real solve_start_time = _current_time;
    for( ; _current_time < solve_start_time + simulation_coupled_time_step; _current_time += simulation_time_step )
    {
        //Delayed Source Term            
        Real delayed_contributions = 0;
        
        //For Each Delayed Group get 
        // dCdt : rate of change
        // delayed_contributions : neutrons that decay out of the group and into the power
        for(int index=0; index < _delayed_precursors.size(); index++ )
        {
            Real lambda_k = _delayed_neutron_set._beta_time_constants[ index ];
            Real precursor_delayed_contribution = _delayed_neutron_set._delayed_neutrons_per_fission_groups[ index ]; 
            Real dCkdt = -lambda_k *  _delayed_precursors[index] + precursor_delayed_contribution * _current_power / _power_per_fission;
            delayed_contributions += lambda_k *  _delayed_precursors[index]; 
            _delayed_precursors[index] += dCkdt * simulation_time_step;
        }
        
        //For every 100000 iterations capture the state of the differential equations for later plotting
        if(iterations % 1000000 == 0)
        {
            std::pair<Real,Real> power = { _current_time, _current_power };
            power_record.push_back(power);
            std::pair<Real,std::vector<Real>> delayed = { _current_time, _delayed_precursors };
            delayed_record.push_back(delayed);
        }
        //Ohhh variable time steps!!!!!!
        
        //Calculate the current power
        Real dPdt = ( reactivity - total_beta )/neutron_generation_time * _current_power + _power_per_fission/(neutron_generation_time * total_neutrons_per_fission)*delayed_contributions;
        _current_power += simulation_time_step * dPdt;
        
        iterations++;
    }
    
    return _current_power;
}
