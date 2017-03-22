/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   infiniteCompositeWorth.h
 * Author: chris
 *
 * Created on March 15, 2017, 4:42 AM
 */

#ifndef INFINITECOMPOSITEWORTH_H
#define INFINITECOMPOSITEWORTH_H

#include "WorthStudy.h"
#include "CompositeMicroCell.h"

class InfiniteCompositeWorth : public WorthStudy
{
public:
    
    CompositeMicroCell* _thermal_solver;
    InfiniteCompositeReactor* _reactor;
    
    InfiniteCompositeWorth(InfiniteCompositeReactor* reactor);
    virtual ~InfiniteCompositeWorth();
    virtual void startStudy(const std::string &output_file_name = "worth.csv");
    virtual void startStudy(const bool &vary_fuel_temperature,const bool &vary_matrix_temperature,const std::string &output_file_name);   
    virtual void setThermalSolver(CompositeMicroCell* thermal_solver);
    
    
    
private:

};

#endif /* INFINITECOMPOSITEWORTH_H */

