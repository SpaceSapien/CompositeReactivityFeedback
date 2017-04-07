/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CylindricalMesh.cpp
 * Author: chris
 * 
 * Created on April 4, 2017, 7:41 PM
 */

#include <vector>
#include "CylindricalMesh.h"
#include <cmath>

CylindricalMesh::CylindricalMesh( const Real &coolant_channel_radius, const Real &outer_radius, const int &number_nodes, const int &number_cells)
{
    
    Dimension radial_length = (outer_radius - coolant_channel_radius) / number_nodes; 
    
    
    for(std::size_t current_node = 1; current_node <= number_nodes; ++current_node )
    {
        
        Dimension current_location = coolant_channel_radius + (current_node - 0.5 ) * radial_length;
        Materials current_material = Materials::Homogenous;
        this->pushBack(current_location, current_material);
            
        Dimension cell_inner_radius = coolant_channel_radius + (current_node - 1 ) * radial_length;
        Dimension cell_outer_radius = coolant_channel_radius + (current_node ) * radial_length;
        
        _inner_radius.push_back(cell_inner_radius);
        _outer_radius.push_back(cell_outer_radius);
            
        Real outer_volume = M_PI * std::pow(cell_outer_radius,2);
        Real inner_volume = M_PI * std::pow(cell_inner_radius,2);
        
        Real volume = outer_volume - inner_volume;
        _volume.push_back(volume);
            
        Real inner_surface_area = 2.0 * M_PI * cell_inner_radius;
        Real outer_surface_area = 2.0 * M_PI * cell_outer_radius;
        
        _inner_surface.push_back(inner_surface_area);
        _outer_surface.push_back(outer_surface_area);            
        
        // Only one zone homogenized
        _zone.push_back(1);
        
        Real range = outer_radius - coolant_channel_radius ;
        Real fraction_range = ( current_location - coolant_channel_radius) / range ;
        int current_cell = std::floor( number_cells * fraction_range) + 1;
        
        _cell_in_zone.push_back(current_cell);       
        
    }   
   
}