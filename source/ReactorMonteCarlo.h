/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactorMonteCarlo.h
 * Author: chris
 *
 * Created on December 17, 2015, 7:09 PM
 */

#ifndef REACTORMONTECARLO_H
#define REACTORMONTECARLO_H
#include <fstream>
#include <iostream>
#include <string>
#include "MicroCell.h"
#include "InfiniteCompositeReactor.h"

class InfiniteCompositeReactor;

class ReactorMonteCarlo 
{

public:
    
    /**
     *
     */
    Real _virtual_k_eff_multiplier;
    Real _current_k_eff;
    Real _current_k_eff_sigma;
    Real _current_prompt_neutron_lifetime;
    Real _current_prompt_neutron_lifetime_sigma;
    
    std::string _run_directory;
    std::vector<std::pair<FissionableIsotope,Real> > _fission_tally_listing;
    
    ReactorMonteCarlo();
    ReactorMonteCarlo(InfiniteCompositeReactor* reactor, const Real &starting_k_effective,const std::string &run_directory);
    void createMCNPOutputFile(const std::string &file_name);
    void updateAdjustedCriticalityParameters();
    void getRawCriticalityParameters( Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma);   
    void readOutputFile(const std::string &file_name, Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma);
    
    std::string getMaterialCards();
    std::string getCellCards();
    std::string getSurfaceCards();
    
    
private:
    
    InfiniteCompositeReactor* _reactor;

};

#endif /* REACTORMONTECARLO_H */

