/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BetaSimulationResults.cpp
 * Author: chris
 * 
 * Created on December 26, 2016, 2:41 PM
 */
#include <cmath>
#include "BetaSimulationResults.h"


BetaSimulationResults::BetaSimulationResults(const SimulationResults &with_delayed_neutrons, const SimulationResults &without_delayed_neutrons) 
{
    _with_delayed_neutrons = SimulationResults(with_delayed_neutrons );
    _without_delayed_neutrons = SimulationResults(without_delayed_neutrons);
    
    _k_eff = _with_delayed_neutrons._k_eff;
    _k_eff_sigma = _with_delayed_neutrons._k_eff_sigma;
    _prompt_neutron_lifetime = _with_delayed_neutrons._prompt_neutron_lifetime;
    _prompt_neutron_lifetime = _with_delayed_neutrons._prompt_neutron_lifetime_sigma;
    
    _beta_sigma = getBetaEffSigma(_k_eff,_k_eff_sigma, _without_delayed_neutrons._k_eff, _without_delayed_neutrons._k_eff_sigma );
    _beta = ( _k_eff - _without_delayed_neutrons._k_eff ) / _k_eff;
    _elapsed_time = with_delayed_neutrons._elapsed_time + without_delayed_neutrons._elapsed_time;
}

/**
 * Helper function to calculate the beta sigma
 * @param k_eff
 * @param k_eff_sigma
 * @param nd_k_eff
 * @param nd_k_eff_sigma
 * @return 
 */
Real BetaSimulationResults::getBetaEffSigma(const Real &k_eff,const Real &k_eff_sigma,const Real &nd_k_eff,const Real &nd_k_eff_sigma)
{
    Real current_beta_eff = ( k_eff - nd_k_eff ) / k_eff;
    
    //sqrt of the sum of the squares of sigmas
    Real difference_uncertainty = sqrt( k_eff_sigma * k_eff_sigma + nd_k_eff_sigma * nd_k_eff_sigma );
    
    //Find the uncertainty fractions
    Real numerator_uncertainty_fraction = difference_uncertainty / ( k_eff - nd_k_eff );
    Real denominator_uncertainty_fraction = k_eff_sigma / k_eff;
    
    //uncertainty is the sum of the squares of the fractions
    return current_beta_eff * sqrt( numerator_uncertainty_fraction * numerator_uncertainty_fraction  + denominator_uncertainty_fraction * denominator_uncertainty_fraction );
}
