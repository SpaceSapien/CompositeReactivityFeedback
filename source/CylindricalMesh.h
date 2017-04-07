/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CylindricalMesh.h
 * Author: chris
 *
 * Created on April 4, 2017, 7:41 PM
 */

#ifndef CYLINDRICALMESH_H
#define CYLINDRICALMESH_H

#include <vector>
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "Mesh.h"
#include "MicroGeometry.h"

class CylindricalMesh : public Mesh
{
public :
    
    //Homogenous case
    CylindricalMesh( const Real &coolant_channel_radius, const Real &outer_radius, const int &number_nodes, const int &number_cells_per_zone);
    
};



#endif /* CYLINDRICALMESH_H */

