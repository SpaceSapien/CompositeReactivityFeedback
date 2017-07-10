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
#include "ReactivityInsertion.h"
#include "InfiniteCompositeReactor.h"
#include "MaterialLibrary.h"
#include <iostream>

ReactorKinetics::ReactorKinetics() { }


/**
 * Constructor to get the initial conditions for the reactor kinetics equations
 * @param initial_power
 * @param state
 */
ReactorKinetics::ReactorKinetics(Reactor* reactor, const Real &initial_power,const DelayedPrecursorInitialState &state, const Real &starting_beta_eff) 
{
    _reactor = reactor;
    _initial_power = initial_power;
    _current_power = _initial_power;
    _current_time  = 0;    
    _kinetics_time_step = _reactor->_input_file_reader->getInputFileParameter("Kinetics Time Iteration", static_cast<Real>(10e-9) );
    _ortensi = _reactor->_input_file_reader->getInputFileParameter("Ortensi", false);
    
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
            
            //The calculated b_eff is going to be a little different than the real b_eff so update the 
            //delayed yields from the 6 groups but use them proportionally   
            Real natural_b_eff = _delayed_neutron_set._delayed_neutrons_per_fission / _delayed_neutron_set._neutrons_per_fission;
            
            Real delayed_neutron_group_proportionality_constant;
            
            if(_ortensi)
            {
                delayed_neutron_group_proportionality_constant = 1;
            }
            else
            {
                delayed_neutron_group_proportionality_constant = starting_beta_eff / natural_b_eff;
            }
            
            //initialize the delated precoursors
            for(int index = 0; index < _delayed_precursors.size(); index++)
            {
                _delayed_neutron_set._delayed_neutrons_per_fission_groups[index] *=  delayed_neutron_group_proportionality_constant;
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


Real ReactorKinetics::getFixedEigenValue(const Real &current_time)
{
    Real k_effective;
    
    if(_ortensi)
    {
        std::vector<Real> time   =  {0, 0.1,   0.25,  0.45,  0.55,  4,  5 };
        std::vector<Real> keigen =  {1, 1.025, 1.018, 1.003,  1.002, .0998 , .997 };

        std::pair<Real,Real> pair = MaterialLibrary::interpolateDataAndTemperatureArraysAndDerivative(time,keigen,current_time);
        k_effective = pair.first;
    }
    else
    {
         k_effective = _reactor->_reactivity_insertion_model->getCurrentKeff(current_time);   
    }
    return k_effective;
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
    const Real &simulation_coupled_time_step
)
{
    
    if(!_ortensi)
    {
    
        Real beta_effective = _reactor->_monte_carlo_model->_current_beta_eff;
        Real total_neutrons_per_fission = _delayed_neutron_set._neutrons_per_fission;    
        Real simulation_time_step = 10e-9;
        //Loop over a simulation coupled time step
        Real solve_start_time = _current_time;
        for( ; _current_time < solve_start_time + simulation_coupled_time_step; _current_time += simulation_time_step )
        {

            Real k_effective = _reactor->_reactivity_insertion_model->getCurrentKeff(_current_time);       
            Real neutron_generation_time = _reactor->_monte_carlo_model->_current_prompt_neutron_lifetime / k_effective;
            Real reactivity = ( k_effective - 1 )/k_effective;

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
    else
    {
       
        std::vector<Real> time =    { 0,0.0489477906245044,0.0802110967178105,0.129170896123136,0.184255173746934,0.221062086961877,0.251768539519987,0.294699930953403,0.331518852949167,0.343863879632678,0.399020209941399,0.448016035689186,0.484726878657565,0.533566590254685,0.600743710164763,0.637370491667397,0.692310663921350,0.771700713925913,0.875528632900267,0.985493038873915,1.13819669578785,1.29093637904425,1.46202547939442,1.82866556662621,2.18921720198198,2.53146745536726,2.87372971753335,3.21599197969945,3.57662767652096,4.02281392790798,4.84797328442947, 100};
        std::vector<Real> power_multiple =  { 1,1.30000000000000,1.60000000000000,3.45776031434188,5.34381139489199,8.80157170923387,12.5736738703340,16.3457760314342,20.1178781925344,23.2612966601179,27.0333988212181,29.5481335952849,30.4911591355600,28.9194499017682,27.3477406679765,26.0903732809431,24.2043222003930,22.3182711198429,20.1178781925344,18.5461689587427,15.7170923379175,13.8310412573674,12.2593320235757,9.43025540275056,7.22986247544208,5.97249508840868,5.02946954813364,4.08644400785859,4.08644400785859,3.45776031434188,2.82907662082519,  1};

        std::pair<Real,Real> pair = MaterialLibrary::interpolateDataAndTemperatureArraysAndDerivative(time,power_multiple,_current_time);
        Real base_power = pair.first;
         _current_time += simulation_coupled_time_step; 
        return base_power * this->_initial_power;
    }
}

