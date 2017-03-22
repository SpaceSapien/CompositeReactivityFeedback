/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SphericalMesh.h
 * Author: chris
 *
 * Created on November 1, 2016, 2:57 PM
 */

#ifndef SPHERICAL_MESH_H
#define SPHERICAL_MESH_H

#include <vector>
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "Mesh.h"
#include "MicroGeometry.h"

class SphericalMesh : public Mesh
{
public :
    
    SphericalMesh(MaterialLibrary::MicroGeometry* geometry,const int &minumum_nodes_per_cell,const int &total_cells, const int &cells_per_zone);
    
};


#endif /* SPHERICAL_MESH_H */

