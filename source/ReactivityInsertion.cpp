/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactivityInsertion.cpp
 * Author: chris
 * 
 * Created on January 31, 2017, 3:35 PM
 */

#include "ReactivityInsertion.h"

const std::string ReactivityInsertion::INSTANTANEOUS_REACTIVITY_INSERTION = "Instantaneous";
const std::string ReactivityInsertion::RAMP_REACTIVITY_INSERTION = "Ramp";



ReactivityInsertion::ReactivityInsertion(Reactor* reactor) 
{
    _reactor = reactor;    
    _initial_k_eff =  this->_reactor->_input_file_reader->getInputFileParameter("Starting K-eff",static_cast<Real>(1.01) );   
    _reactivity_insertion_function = this->_reactor->_input_file_reader->getInputFileParameter("Reactivity Insertion Method", ReactivityInsertion::INSTANTANEOUS_REACTIVITY_INSERTION);
    _ending_virtual_k_eff_multiplier = _initial_k_eff / _reactor->_monte_carlo_model->_current_k_eff;
    
    if( _reactivity_insertion_function == ReactivityInsertion::RAMP_REACTIVITY_INSERTION)
    {
        //If not specified make the ramp over 1 ms
        _reactivity_addition_timing = this->_reactor->_input_file_reader->getInputFileParameter("Ramp Insertion Timing", 0.001);
        _starting_virtual_k_eff_multiplier = 1.0 / _reactor->_monte_carlo_model->_current_k_eff;
        _ramp_finished = false;
    }
    else
    {
        _reactivity_addition_timing = -1;
        _starting_virtual_k_eff_multiplier = _ending_virtual_k_eff_multiplier;     
        _ramp_finished = true;
    }
}

//A function that is called to break the inner loop and record a timestep right as the ramp ends
bool ReactivityInsertion::rampNeedsReactivityMonteCarloUpdate(const Real &last_update_time,const Real &time_since_last_update )
{
    Real current_time = last_update_time + time_since_last_update;
    
    //if this is a ramp
    if( _reactivity_insertion_function == ReactivityInsertion::RAMP_REACTIVITY_INSERTION )
    {
        //the ramp finished flag hasn't been set
        if(! _ramp_finished)
        {
            //and the current time is over the ramp addition time
            if( current_time >= _reactivity_addition_timing )
            {
                _ramp_finished = true;
                return true;
            }
        }
    }
    return false;
}

bool ReactivityInsertion::rampNeedsReactivityLogUpdate(const Real &last_update_time,const Real &time_since_last_update )
{
    Real current_time = last_update_time + time_since_last_update;
    
    //if this is a ramp
    if( _reactivity_insertion_function == ReactivityInsertion::RAMP_REACTIVITY_INSERTION )
    {
        //the ramp finished flag hasn't been set
        if(! _ramp_finished)
        {
            //There needs to be at leat 25 piecewise time periods where the reactivity needs to be shown in the log 
            if( time_since_last_update > _reactivity_addition_timing/25.0)
            {
                return true;
            }
        }
    }
    return false;
}

Real ReactivityInsertion::getCurrentVirtualKeffMultiplier(const Real &current_time)
{
    if( _reactivity_addition_timing >= 0.0 && current_time < _reactivity_addition_timing)
    {
        Real delta_multiplier = _ending_virtual_k_eff_multiplier - _starting_virtual_k_eff_multiplier;        
        Real fraction_time_in = current_time/_reactivity_addition_timing;
        Real multiplier = _starting_virtual_k_eff_multiplier + fraction_time_in * delta_multiplier;
        return multiplier;
    }
    else
    {
        return _ending_virtual_k_eff_multiplier;
    }  
}

Real ReactivityInsertion::getCurrentKeff(const Real &current_time)
{
    Real multiplier = this->getCurrentVirtualKeffMultiplier(current_time);
    Real current_k_eff = _reactor->_monte_carlo_model->_current_k_eff;
    return current_k_eff * multiplier;
}

Real ReactivityInsertion::getCurrentKeffSigma(const Real &current_time)
{
    Real multiplier = this->getCurrentVirtualKeffMultiplier(current_time);
    Real current_k_eff_sigma = _reactor->_monte_carlo_model->_current_k_eff_sigma;
    return current_k_eff_sigma * multiplier;
}

