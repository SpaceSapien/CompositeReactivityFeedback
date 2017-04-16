/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Reactor.h
 * Author: chris
 *
 * Created on March 12, 2017, 10:34 PM
 */

#ifndef REACTOR_H
#define REACTOR_H

#include <string>
#include <map>
#include <ctime>
#include <vector>
#include "EnumsAndFunctions.h"
#include "InputFileParser.h"
#include "MicroSolution.h"
#include "ReactorMonteCarlo.h"
#include "MaterialLibrary.h"
//#include "ReactivityInsertion.h"

/*
#include "MicroCell.h"
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "InputDataFunctions.h"

#include "PythonPlot.h"
#include "WorthStudy.h"*/

class ReactorMonteCarlo;
class MicroCell;
class MicroGeometry;
class ReactivityInsertion;
class ReactorKinetics;
//class MicroCell;

using namespace MaterialLibrary;

class Reactor 
{
    public:
    
    //For determining a time based or temperature based MC timestep
    enum MoneCarloRecalculation
    {
        Time,
        Temperature
    };  

    MoneCarloRecalculation static getRecalculationType(const std::string &string_type);
    MoneCarloRecalculation _monte_carlo_reclaculation_type;
    Real _maximum_allowed_temperature_difference;
    Real _maximum_allowed_average_temperature_difference;


    //Simulation start time in real time
    const static std::time_t _simulation_start_time;
    static bool _otf_sab;
    
     //The geometry object
    MaterialLibrary::MicroGeometry* _micro_sphere_geometry;
    //Create the thermal heat transfer object 
    MicroCell* _thermal_solver;
    //Create the kinetics model object
    ReactorKinetics* _kinetics_model;
    //Gather the monte carlo object
    ReactorMonteCarlo* _monte_carlo_model;
    //Reactivity Insertion Control for Virtual Insertions
    ReactivityInsertion* _reactivity_insertion_model;    
    //Object to read the input file
    InputFileParser* _input_file_reader;
    //Output file data
    std::string _results_directory;
    std::string _run_name;
    std::string _data_file;            //Data file


    //Time stepping parameters
    Real _transient_time;              //Time since the start of the transient
    Real _monte_carlo_time_iteration;  //How often to calculate keff and the prompt neutron lifetime
    Real _kinetics_thermal_sync_time_step;     //How often to couple the kinetics and heat transfer routines    
    Real _end_time;                    //How many seconds should the simulation last?
    Real _power_and_delayed_neutron_record_time_step; //How often to record the power and delated neutron data
    Real _inner_time_step;
    int _monte_carlo_number_iterations;
    

    //Data Storage for various data
    std::vector<MicroSolution> _plot_solutions;    
    std::vector<std::pair<Real,Real>> _power_record;
    std::vector< std::pair<Real,std::vector<Real>> > _delayed_record;
    std::vector<std::tuple<Real,Real,Real>> _k_eff_record;
    std::vector<std::tuple<Real,Real,Real>> _beta_eff_record;
    std::vector<std::tuple<Real,Real,Real>> _prompt_life_time_record;
    std::vector<std::tuple<Real,Real,Real>> _reactivity_pcm_record;
    std::vector<std::tuple<Real,Real,Real>> _reactivity_cents_record;
    
    //Constructors
    Reactor(const std::string &input_file_name);
    //Reactor(const std::string old_results_folder, Real new_end_time);    
    virtual ~Reactor();
    
    //Implemented but override-able 
    virtual void simulateTransient();
    virtual bool significantTemperatureDifference(MicroSolution* comparison);
    virtual void initializeReactorProblem();
    virtual void monteCarloTimeStepSimulationDataProcessing();
    
    
    //Need to be implemented in base class
    virtual void worthStudy() = 0;
    //virtual void timeIterationInnerLoop() = 0;
    //virtual void temperatureIterationInnerLoop() = 0;
    
    virtual void temperatureIterationInnerLoop();
    virtual void postSimulationProcessing();
    virtual void timeIterationInnerLoop();
    
    
    
    virtual void saveCurrentData
    (
        const Real &time, 
        const Real &power, 
        const Real &k_eff, 
        const Real &k_eff_sigma, 
        const Real &neutron_lifetime, 
        const Real &neutron_lifetime_sigma, 
        const Real &beta_eff, 
        const Real &beta_eff_sigma, 
        const Real &hot_temperature, 
        const Real &gamma, 
        const Real &integrated_power_out, 
        const Real &integrated_power, 
        const Real &current_power_out
    );
    virtual void createOutputFile();
    
    
    virtual void setThermalSolver(MicroCell* solver);
    virtual void setMicroSphereGeometry(MaterialLibrary::MicroGeometry* geometry);
    virtual void setKineticsModel(ReactorKinetics* kinetics_model);
    virtual void setMoteCarloModel(ReactorMonteCarlo* monte_carlo);
    virtual void setReactivityInsertionModel(ReactivityInsertion* insertion_model);
    virtual void setInputFileParser(InputFileParser* file_reader);
    
      
    
    
    private:

};

#endif /* REACTOR_H */

