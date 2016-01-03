/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ExplicitSolverSettings.h
 * Author: chris
 *
 * Created on November 20, 2015, 5:49 PM
 */

#ifndef EXPLICITSOLVERSETTINGS_H
#define EXPLICITSOLVERSETTINGS_H
#include "EnumsAndFunctions.h"
#include <vector>

class ExplicitSolverSettings 
{
public:

    enum SolverOrder
    {
        SECOND=1,  //the numbers denot the grid points to the left and right
        FOURTH=2   //that are used for these solver orders
    };
    
    ExplicitSolverSettings();
    ExplicitSolverSettings( const SolverOrder &solver_order, const Dimension &outer_radius); 
    void printSolverSettings();
    
    Real _time_step;
    Real _start_time;
    Real _max_radius;
    int _number_elements;
    Real _element_size;
    int _number_points;
    SolverOrder _solver_order;
    std::vector<Real> _radial_mesh;
    
private:

};

#endif /* EXPLICITSOLVERSETTINGS_H */

