/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SimulationResults.h
 * Author: chris
 *
 * Created on March 1, 2016, 1:44 PM
 */

#ifndef SIMULATIONRESULTS_H
#define SIMULATIONRESULTS_H
#include <string>
#include <vector>
#include "EnumsAndFunctions.h"


class TallyResults 
{
    
public:
    
    void readTallyFile(const std::string &output_file, const std::string &directory);
    TallyResults(const std::string &output_file, const std::string &directory );
    TallyResults();
    
private:

};


/**
 * This class represents the data that can be represented from an output file
 * 
 * 
 */
class SimulationResults 
{
public:
    
    Real _k_eff;
    Real _k_eff_sigma;
    Real _prompt_neutron_lifetime;
    Real _prompt_neutron_lifetime_sigma;
    
    void readOutputFile(const std::string &output_file, const std::string &directory);
    SimulationResults(const std::string &output_file, const std::string &directory );
    SimulationResults();
    SimulationResults(const SimulationResults &results);
    
private:

};

class BetaSimulationResults : SimulationResults
{
public: 
    
    Real _beta;
    Real _beta_sigma;
    
    SimulationResults _with_delayed_neutrons;
    SimulationResults _without_delayed_neutrons;
    
    BetaSimulationResults(const SimulationResults &with_delayed_neutrons, const SimulationResults &without_delayed_neutrons);
    static Real getBetaEffSigma(const Real &k_eff,const Real &k_eff_sigma,const Real &nd_k_eff,const Real &nd_k_eff_sigma);
    
};


#endif /* SIMULATIONRESULTS_H */

