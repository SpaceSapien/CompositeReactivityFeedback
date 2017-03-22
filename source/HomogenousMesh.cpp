/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMesh.cpp
 * Author: chris
 * 
 * Created on March 16, 2017, 3:58 PM
 */

#include "HomogenousMesh.h"
#include "MaterialLibrary.h"



HomogenousMesh::HomogenousMesh(MaterialLibrary::MicroGeometry* geometry) 
{
    _inner_radius.push_back(0);
    _outer_radius.push_back(1);   
    _volume.push_back(1);
    _inner_surface.push_back(0);
    _outer_surface.push_back(6);            
    _zone.push_back(1);
    _cell_in_zone.push_back(1);
    _position.push_back(0.0);
    _materials.push_back(Materials::Homogenous);
}


