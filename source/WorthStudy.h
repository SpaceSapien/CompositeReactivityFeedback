/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   WorthStudy.h
 * Author: chris
 *
 * Created on December 17, 2016, 7:47 PM
 */

#ifndef WORTHSTUDY_H
#define WORTHSTUDY_H
#include <string>
#include "MicroCell.h"
#include "MicroGeometry.h"
#include "ReactorMonteCarlo.h"

class WorthStudy 
{
public:
    
    //The geometry object
    MicroGeometry* _micro_sphere_geometry;
    //Create the thermal heat transfer object 
    ReactorMonteCarlo* _monte_carlo_model;
    //Read the input file
    InfiniteCompositeReactor* _reactor;
    //The thermal data and mesh
    MicroCell* _thermal_solver;
    std::string _output_file;
    
    
    WorthStudy(InfiniteCompositeReactor* reactor);
    void startStudy();
    void log(const SimulationResults &results);
    void createOutputFile();
    
private:

};

#endif /* WORTHSTUDY_H */

