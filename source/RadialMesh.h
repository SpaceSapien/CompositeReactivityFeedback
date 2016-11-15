/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   RadialMesh.h
 * Author: chris
 *
 * Created on November 1, 2016, 2:57 PM
 */

#include <vector>
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"

#ifndef RADIALMESH_H
#define RADIALMESH_H

class RadialMesh
{
public :
    
    std::vector<Dimension> _position;
    std::vector<Materials> _materials;
    std::vector<Real> _volume;
    std::vector<Real> _inner_radius;
    std::vector<Real> _outer_radius;
    std::vector<Real> _inner_surface;
    std::vector<Real> _outer_surface;
    std::vector<int> _zone;
    std::vector<int> _cell_in_zone;
    
    RadialMesh(MicroGeometry* geometry,const int &minumum_nodes_per_cell,const int &total_cells, const int &cells_per_zone);
    
    std::size_t numberOfNodes() const;        
    void pushBack(const Dimension &dimension, const Materials &material);
    Dimension getNodeLocation(const int &node) const;
    Materials getMaterial(const int &node) const;
    
};


#endif /* RADIALMESH_H */

