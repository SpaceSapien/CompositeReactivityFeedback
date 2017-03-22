/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMonteCarlo.h
 * Author: chris
 *
 * Created on March 17, 2017, 12:55 PM
 */

#ifndef HOMOGENOUSMONTECARLO_H
#define HOMOGENOUSMONTECARLO_H

#include "ReactorMonteCarlo.h"
#include "InfiniteHomogenousReactor.h"

class HomogenousMonteCarlo : public ReactorMonteCarlo
{
public:
    
    HomogenousMonteCarlo(InfiniteHomogenousReactor* reactor, const std::string &run_directory);    
    virtual ~HomogenousMonteCarlo();
    
    virtual std::string getMaterialCards();
    virtual std::string getCellCards();
    virtual std::string getSurfaceCards();
    
    
private:
    
    InfiniteHomogenousReactor* _reactor;

};

#endif /* HOMOGENOUSMONTECARLO_H */

