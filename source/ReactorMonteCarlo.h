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
#include "SimulationResults.h"
#include "TallyGroup.h"
#include "BetaSimulationResults.h"

class InfiniteCompositeReactor;

class ReactorMonteCarlo 
{

public:
    
    
    Real _current_k_eff;
    Real _current_k_eff_sigma;
    Real _current_prompt_neutron_lifetime;
    Real _current_prompt_neutron_lifetime_sigma;
    Real _current_beta_eff;
    Real _current_beta_eff_sigma;
    long _current_number_particles;
    
    Real _virtual_k_eff_multiplier;
    Real _starting_k_eff;
    int _number_cpus;
    int _cells_per_zone;
    int _number_zones;
    
    long _beta_eff_number_particles;
    long _k_eff_number_particles;
    long _k_eff_number_cycles;
    long _beta_eff_number_cycles;
    long _particles_per_cycle;
    
    
    int _calulate_beta_interval;
    int _number_of_keff_calculations;
    std::time_t _current_mc_exection_elapsed_time;
    
    
    
    
    bool _tally_cells;
    int _tally_energy_bins;
    
    std::vector<TallyGroup*> _tally_groups;
    
    std::string _run_directory;
    
   
    ~ReactorMonteCarlo();
    ReactorMonteCarlo(InfiniteCompositeReactor* reactor, const Real &starting_k_effective,const std::string &run_directory);
    void createMCNPOutputFile(const std::string &run_title, const std::string &file_name,const int &particles_per_cycle,const int &number_cycles, const bool &delated_neutrons = true);
    void updateAdjustedCriticalityParameters();
    
    void updateCurrentValuesFromResults(const BetaSimulationResults &results);
    void updateCurrentValuesFromResults(const SimulationResults &results);
    
    BetaSimulationResults getRawKeffAndBetaEff();
    SimulationResults getRawKeff();    
    SimulationResults getRawCriticalityParameters(const std::string &file_root, const int &particles_per_cycle, const int &number_cycles, const bool &delayed_neutrons);
    
    void readOutputFile(const std::string &file_name, Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma);
    
    std::string getMaterialCards();
    std::string getCellCards();
    std::string getSurfaceCards();
    std::string getTallyCards();
    std::string getSingleCellCard(const Materials &material, const int &current_zone, int &cell_number );
    
    void createTallyOutputFile(std::string file_name = "tally-data.csv");
    void outputTalliesToFile(TallyGroup* tally_group, std::string file_name = "tally-data.csv");
    
    std::vector< std::vector<Real> > getZoneCellRelativePowerDensity();
    
private:
    
    //Helper function to calculate uncertainties
    //Pointer the the reactor
    InfiniteCompositeReactor* _reactor;

};

#endif /* REACTORMONTECARLO_H */

