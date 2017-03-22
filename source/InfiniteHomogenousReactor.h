/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InfiniteHomogenousReactor.h
 * Author: chris
 *
 * Created on March 13, 2017, 11:54 AM
 */

#ifndef INFINITEHOMOGENOUSREACTOR_H
#define INFINITEHOMOGENOUSREACTOR_H
#include <string>
#include "Reactor.h"
#include "HomogenousMicroCell.h"


class HomogenousMicroCell;

class InfiniteHomogenousReactor : public Reactor
{
public:
    
    
    //Create the thermal heat transfer object 
    HomogenousMicroCell* _thermal_solver;
    
    InfiniteHomogenousReactor(const std::string &name);
    
    //virtual void timeIterationInnerLoop();
    //virtual void temperatureIterationInnerLoop();
    virtual void initializeReactorProblem();
    virtual void worthStudy();
    //virtual void postSimulationProcessing();    
    //virtual void monteCarloTimeStepSimulationProcessing();
    
    virtual void setThermalSolver(HomogenousMicroCell* solver);
    
    virtual ~InfiniteHomogenousReactor();
private:

};

#endif /* INFINITEHOMOGENOUSREACTOR_H */

