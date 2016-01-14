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
#include <math.h>
#include <iomanip>
#include <ctime>
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "ExplicitSolverSettings.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"



class InfiniteCompositeReactor
{
    public:    
    
    //The geometry object
    MicroGeometry* _micro_sphere_geometry;
    //Create the thermal heat transfer object 
    MicroCell* _thermal_solver;
    //Create the kinetics model object
    ReactorKinetics* _kinetics_model;
    //Gather the monte carlo object
    ReactorMonteCarlo* _monte_carlo_model;
    
    //Data Storage for various data
    std::vector<MicroSolution> _plot_solutions;    
    std::vector<std::pair<Real,Real>> _power_record;
    std::vector< std::pair<Real,std::vector<Real>> > _delayed_record;
    std::vector<std::tuple<Real,Real,Real>> _k_eff_record;
    std::vector<std::tuple<Real,Real,Real>> _prompt_life_time_record;
    
    
    //Time stepping parameters
    Real _monte_carlo_time_iteration;  //How often to calculate keff and the prompt neutron lifetime
    Real _kinetics_time_iteration;     //How often to couple the kinetics and heat transfer routines    
    Real _end_time;                    //How many seconds should the simulation last?
    std::string _data_file;            //Data file
    
    InfiniteCompositeReactor();
    virtual ~InfiniteCompositeReactor();
    std::string getSaveDirectory();
    std::string _results_directory;
    
    void simulate();
    void initializeInifiniteCompositeReactorProblem();
    void plotDelayedPrecursors();
    void saveCurrentData(const Real &time, const Real &power, const Real &k_eff, const Real &k_eff_sigma, const Real &neutron_lifetime, const Real &neutron_lifetime_sigma);
    void createOutputFile();
private:
    
    
};

#endif 

