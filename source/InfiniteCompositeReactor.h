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
#include "Reactor.h"
#include "InfiniteCompositeMonteCarlo.h"
/*#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "WorthStudy.h"
#include "ReactivityInsertion.h"*/

class CompositeMicroCell;
class InfiniteCompositeMonteCarlo;

class InfiniteCompositeReactor : public Reactor
{
public:    
     
    //Create the thermal heat transfer object 
    CompositeMicroCell* _thermal_solver;
    InfiniteCompositeMonteCarlo* _monte_carlo_model;
    
    InfiniteCompositeReactor(const std::string &input_file = "");
    InfiniteCompositeReactor(const std::string resume_file, Real new_end_time);
    virtual ~InfiniteCompositeReactor();    
    
    //virtual void timeIterationInnerLoop();
    //virtual void temperatureIterationInnerLoop();
    virtual void initializeReactorProblem();
    virtual void worthStudy();
    
    void solveForSteadyStatePowerDistribution(const std::vector<Real> &homogenous_power_density, const Real &initial_power_density);
    
    virtual void setThermalSolver(CompositeMicroCell* solver);
    virtual void setMoteCarloModel(InfiniteCompositeMonteCarlo* mc_model);
    
    
protected:
    
    
};

#endif 

