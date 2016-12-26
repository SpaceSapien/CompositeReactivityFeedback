/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BetaSimulationResults.h
 * Author: chris
 *
 * Created on December 26, 2016, 2:41 PM
 */

#ifndef BETASIMULATIONRESULTS_H
#define BETASIMULATIONRESULTS_H
#include "SimulationResults.h"

class BetaSimulationResults : public SimulationResults
{
public: 
    
    Real _beta;
    Real _beta_sigma;
    
    SimulationResults _with_delayed_neutrons;
    SimulationResults _without_delayed_neutrons;
    
    BetaSimulationResults(const SimulationResults &with_delayed_neutrons, const SimulationResults &without_delayed_neutrons);
    static Real getBetaEffSigma(const Real &k_eff,const Real &k_eff_sigma,const Real &nd_k_eff,const Real &nd_k_eff_sigma);
    
   
};

#endif /* BETASIMULATIONRESULTS_H */

