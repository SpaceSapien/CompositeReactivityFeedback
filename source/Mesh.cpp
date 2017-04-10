/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mesh.cpp
 * Author: chris
 * 
 * Created on March 16, 2017, 3:47 PM
 */

#include "Mesh.h"

void Mesh::pushBack(const Dimension &dimension, const Materials &material)
{
     _materials.push_back( material );
     _position.push_back( dimension );
}

Dimension Mesh::getNodeLocation(const int &node) const
{
    return _position[node];
}

Materials Mesh::getMaterial(const int &node) const
{
    return _materials[node];
}

std::size_t Mesh::numberOfNodes() const
{
    return _position.size();
}

Real Mesh::getTotalVolume() const
{
    std::size_t number_nodes = numberOfNodes();
    Real total_volume = 0;
    
    for(std::size_t node_index = 0; node_index < number_nodes; ++node_index)
    {
        total_volume +=  _volume[node_index];
    }
}

Real Mesh::getVolumeWeightedQuantity(const std::vector<Real> &quantity) const
{
    std::size_t number_nodes = numberOfNodes();
    std::size_t quantity_size = quantity.size();
    
    if(number_nodes != quantity_size)
    {
        throw std::string("Wrong number of nodes");
    }
    
    Real weighted_quantity = 0;
    
    for(std::size_t node_index = 0; node_index < number_nodes; ++node_index)
    {
        weighted_quantity += quantity[node_index] * _volume[node_index];
    }
    
    return weighted_quantity;
}