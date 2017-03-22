/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactivityInsertion.h
 * Author: chris
 *
 * Created on January 31, 2017, 3:35 PM
 */
#include "EnumsAndFunctions.h"
#include "InfiniteCompositeReactor.h"

#ifndef REACTIVITYINSERTION_H
#define REACTIVITYINSERTION_H

//class Reactor;

class ReactivityInsertion 
{
public:
    ReactivityInsertion(Reactor* reactor);
   
    Real getCurrentVirtualKeffMultiplier(const Real &current_time);
    Real getCurrentKeff(const Real &current_time);
    Real getCurrentKeffSigma(const Real &current_time);
    bool rampNeedsReactivityLogUpdate(const Real &last_update_time,const Real &time_since_last_update );
    bool rampNeedsReactivityMonteCarloUpdate(const Real &last_update_time,const Real &time_since_last_update );
    
    
    Reactor* _reactor;
    
    std::string _reactivity_insertion_function;
    Real _initial_k_eff;
    Real _starting_virtual_k_eff_multiplier;
    Real _ending_virtual_k_eff_multiplier;
    Real _reactivity_addition_timing;
    bool _ramp_finished;
    
    static const std::string INSTANTANEOUS_REACTIVITY_INSERTION;
    static const std::string RAMP_REACTIVITY_INSERTION;
    
    
private:

};

#endif /* REACTIVITYINSERTION_H */

