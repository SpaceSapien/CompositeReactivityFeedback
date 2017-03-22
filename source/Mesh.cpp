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