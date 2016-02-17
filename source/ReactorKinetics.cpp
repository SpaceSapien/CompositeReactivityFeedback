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
#include "InfiniteCompositeReactor.h"
#include <iostream>

ReactorKinetics::ReactorKinetics() { }


/**
 * Constructor to get the initial conditions for the reactor kinetics equations
 * @param initial_power
 * @param state
 */
ReactorKinetics::ReactorKinetics(InfiniteCompositeReactor* reactor, const Real &initial_power,const DelayedPrecursorInitialState &state) 
{
    _reactor = reactor;
    _initial_power = initial_power;
    _current_power = _initial_power;
    _current_time  = 0;    
    _kinetics_time_step = _reactor->_input_file_reader->getInputFileParameter("Kinetics Time Iteration", 10e-9);
    
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
    const Real &beta_effective
)
{
    
    //Calculate some parameters
    Real reactivity = ( k_effective - 1 )/k_effective;
    Real total_neutrons_per_fission = _delayed_neutron_set._neutrons_per_fission;
    
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
        
        //Calculate the current power
        Real dPdt = ( reactivity - beta_effective )/neutron_generation_time * _current_power + _power_per_fission/(neutron_generation_time * total_neutrons_per_fission)*delayed_contributions;
        _current_power += simulation_time_step * dPdt;
        
        
    }
    
    return _current_power;
}
