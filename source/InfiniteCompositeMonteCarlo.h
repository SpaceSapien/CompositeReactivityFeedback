/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InfiniteCompositeMonteCarlo.h
 * Author: chris
 *
 * Created on March 17, 2017, 9:08 AM
 */

#ifndef INFINITECOMPOSITEMONTECARLO_H
#define INFINITECOMPOSITEMONTECARLO_H
#include <string>
#include "ReactorMonteCarlo.h"
#include "InfiniteCompositeReactor.h"

class InfiniteCompositeReactor;

class InfiniteCompositeMonteCarlo : public ReactorMonteCarlo
{
public:
    
    
    
    InfiniteCompositeMonteCarlo(InfiniteCompositeReactor* reactor, const std::string &run_dir);    
    virtual ~InfiniteCompositeMonteCarlo();
    
 
    
    virtual std::string getMaterialCards();
    virtual std::string getCellCards();
    virtual std::string getSurfaceCards();
    virtual std::string getSingleCellCard(const Materials &material, const int &current_zone, int &cell_number );
    
private:
    
    InfiniteCompositeReactor* _reactor;

};

#endif /* INFINITECOMPOSITEMONTECARLO_H */

