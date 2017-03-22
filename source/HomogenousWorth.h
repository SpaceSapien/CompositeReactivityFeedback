/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousWorth.h
 * Author: chris
 *
 * Created on March 15, 2017, 7:09 AM
 */

#ifndef HOMOGENOUSWORTH_H
#define HOMOGENOUSWORTH_H
#include "WorthStudy.h"
#include "InfiniteHomogenousReactor.h"
#include "HomogenousMicroCell.h"
#include "HomogenousMonteCarlo.h"

class HomogenousWorth : public WorthStudy
{
public:
    
    InfiniteHomogenousReactor* _reactor;
    HomogenousMicroCell* _thermal_solver;
    HomogenousMonteCarlo* _monte_carlo_model;
    
    HomogenousWorth(InfiniteHomogenousReactor* reactor);
    virtual ~HomogenousWorth();    
    virtual void startStudy(const std::string &output_file_name = "worth.csv");
    
    virtual void setThermalSolver(HomogenousMicroCell* solver);
    virtual void setMonteCarloModel(HomogenousMonteCarlo* monte_carlo_model);
    
    
private:

};

#endif /* HOMOGENOUSWORTH_H */

