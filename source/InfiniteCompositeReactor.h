/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ResultsStorage.h
 * Author: chris
 *
 * Created on December 30, 2015, 3:01 PM
 */

#ifndef INFINITECOMPOSITEREACTOR_H
#define INFINITECOMPOSITEREACTOR_H
#include <tuple>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include <cmath>
#include <iomanip>
#include <ctime>
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InputFileParser.h"

class MicroCell;
class ReactorKinetics;
class ReactorMonteCarlo;


class InfiniteCompositeReactor
{
    public:    
    
        
    enum MoneCarloRecalculation
    {
        Time,
        Temperature
    };  
    
    MoneCarloRecalculation static getRecalculationType(const std::string &string_type);
    MoneCarloRecalculation _monte_carlo_reclaculation_type;
    Real _maximum_allowed_temperature_difference;
    
        
    //The geometry object
    MicroGeometry* _micro_sphere_geometry;
    //Create the thermal heat transfer object 
    MicroCell* _thermal_solver;
    //Create the kinetics model object
    ReactorKinetics* _kinetics_model;
    //Gather the monte carlo object
    ReactorMonteCarlo* _monte_carlo_model;
    //Read the input file
    InputFileParser* _input_file_reader;
    
    //Data Storage for various data
    std::vector<MicroSolution> _plot_solutions;    
    std::vector<std::pair<Real,Real>> _power_record;
    std::vector< std::pair<Real,std::vector<Real>> > _delayed_record;
    std::vector<std::tuple<Real,Real,Real>> _k_eff_record;
    std::vector<std::tuple<Real,Real,Real>> _beta_eff_record;
    std::vector<std::tuple<Real,Real,Real>> _prompt_life_time_record;
    std::vector<std::tuple<Real,Real,Real>> _reactivity_pcm_record;
    std::vector<std::tuple<Real,Real,Real>> _reactivity_cents_record;
    
    //Time stepping parameters
    Real _transient_time;              //Time since the start of the transient
    Real _monte_carlo_time_iteration;  //How often to calculate keff and the prompt neutron lifetime
    Real _kinetics_thermal_sync_time_step;     //How often to couple the kinetics and heat transfer routines    
    Real _end_time;                    //How many seconds should the simulation last?
    Real _power_and_delayed_neutron_record_time_step; //How often to record the power and delated neutron data
    Real _inner_time_step;
    int _monte_carlo_number_iterations;
    std::string _data_file;            //Data file
    
    
    InfiniteCompositeReactor(const std::string &input_file = "");
    virtual ~InfiniteCompositeReactor();
    std::string getSaveDirectory();
    std::string _results_directory;
    std::string _run_name;
    
    
    const static std::time_t _simulation_start_time;
    
    void simulate();
    void timeIterationInnerLoop();
    void temperatureIterationInnerLoop();
    bool significantTemperatureDifference(MicroSolution* comparison);
    
    void initializeInifiniteCompositeReactorProblem();
    void plotDelayedPrecursors();
    void saveCurrentData(const Real &time, const Real &power, const Real &k_eff, const Real &k_eff_sigma, const Real &neutron_lifetime, const Real &neutron_lifetime_sigma, const Real &beta_eff, const Real &beta_eff_sigma, const Real &hot_temperature);
    void createOutputFile();
    void postSimulationProcessing();
    void monteCarloTimeStepSimulationProcessing();
    
    
private:
    
    
};

#endif 

