/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ExplicitSolverSettings.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 5:49 PM
 */
#include <iostream>
#include <vector>
#include "ExplicitSolverSettings.h"


ExplicitSolverSettings::ExplicitSolverSettings(const SolverOrder &solver_order, const Dimension &outer_radius) 
{
    _time_step = 100e-9;
    _start_time = 0;
    _max_radius = outer_radius;
    _number_elements = 100;
    _element_size = _max_radius / _number_elements;
    _number_points = _number_elements + 1;
    _solver_order = solver_order;
    
    _radial_mesh = std::vector<Real>();
    _radial_mesh.reserve( _number_points);
    
    for( int index = 0; index < _number_points; ++index)
    {
        Real r_index = index * _element_size;
        _radial_mesh.push_back( r_index );
    }
    
}

ExplicitSolverSettings::ExplicitSolverSettings()
{}

void ExplicitSolverSettings::printSolverSettings()
{
    std::cout << "Solver Settings @" << this << "\n";
    std::cout << "Time Step " << _time_step << " s \n";
    std::cout << "Max Radius " << _max_radius << " m \n";
    std::cout << "Number Elements " << _number_elements << "  \n";
    std::cout << "Element Size " << _element_size << " m  \n";
    std::cout << "Nodes " << _number_points << "  \n";
    
    
    
    std::cout << "\n";
    
}



