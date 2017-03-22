/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMesh.h
 * Author: chris
 *
 * Created on March 16, 2017, 3:58 PM
 */

#ifndef HOMOGENOUSMESH_H
#define HOMOGENOUSMESH_H
#include "Mesh.h"
#include "MicroGeometry.h"


class HomogenousMesh : public Mesh
{
public:
    HomogenousMesh(MaterialLibrary::MicroGeometry* geometry);

};

#endif /* HOMOGENOUSMESH_H */

