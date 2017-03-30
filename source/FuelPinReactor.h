/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinReactor.h
 * Author: chris
 *
 * Created on March 22, 2017, 10:26 PM
 */

#ifndef FUELPINREACTOR_H
#define FUELPINREACTOR_H

#include "Reactor.h"

class FuelPinReactor : public Reactor
{

public:
    //Create the thermal heat transfer object 
    //HomogenousMicroCell* _thermal_solver;
    
    FuelPinReactor(const std::string &name);
    
    //virtual void timeIterationInnerLoop();
    //virtual void temperatureIterationInnerLoop();
    virtual void initializeReactorProblem();
    virtual void worthStudy();
    
    //virtual void setThermalSolver(HomogenousMicroCell* solver);
    
    virtual ~FuelPinReactor();
    
private:

};

#endif /* FUELPINREACTOR_H */

