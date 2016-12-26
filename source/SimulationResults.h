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
#include <ctime>
#include <string>
#include <vector>
#include "EnumsAndFunctions.h"


/**
 * This class represents the data that can be represented from an output file
 */
class SimulationResults 
{
public:
    
    Real _k_eff;
    Real _k_eff_sigma;
    Real _prompt_neutron_lifetime;
    Real _prompt_neutron_lifetime_sigma;
    std::time_t _elapsed_time;
    
    void readOutputFile(const std::string &output_file, const std::string &directory);
    SimulationResults(const std::string &output_file, const std::string &directory, const std::time_t &elapsed_time );
    SimulationResults(const SimulationResults &results);
    SimulationResults();
    
    
private:

};




#endif /* SIMULATIONRESULTS_H */

