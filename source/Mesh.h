/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mesh.h
 * Author: chris
 *
 * Created on March 16, 2017, 3:47 PM
 */

#ifndef MESH_H
#define MESH_H
#include <vector>
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"

using namespace MaterialLibrary;

class Mesh
{
public:
    
    std::vector<Dimension> _position;
    std::vector<Materials> _materials;
    std::vector<Real> _volume;
    std::vector<Real> _inner_radius;
    std::vector<Real> _outer_radius;
    std::vector<Real> _inner_surface;
    std::vector<Real> _outer_surface;
    std::vector<int> _zone;
    std::vector<int> _cell_in_zone;
    
    
    
    std::size_t numberOfNodes() const;        
    void pushBack(const Dimension &dimension, const Materials &material);
    Dimension getNodeLocation(const int &node) const;
    Materials getMaterial(const int &node) const;
    
    Real getTotalVolume() const;
    Real getVolumeWeightedQuantity(const std::vector<Real> &quantity) const;
    


};

#endif /* MESH_H */

