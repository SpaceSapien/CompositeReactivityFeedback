/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaterialEnum.h
 * Author: chris
 *
 * Created on November 20, 2015, 3:18 PM
 */
#include <string>
#include <vector>
#include <string.h>

#ifndef MATERIALENUM_H
#define MATERIALENUM_H

typedef double Real;
typedef double Dimension;

enum Materials
{
    U,
    UO2,
    UN,
    UC,
    U3Si,
    SiC,
    C,
    Be,
    BeO,
    ZrB2,
    W,
    B4C,
    Mo,
    Nb,
    Zr,
    ZrO2,
    Graphene
};
   
enum FissionableIsotope
{
    U233,
    U235,
    U238,
    Th232,
    Pu239,
    Pu240,
    Pu241,
    Pu242,
    Mixed
};

std::string getMaterialName(const Materials &material);
Materials getMaterialFromName(const std::string &material_name);

Real sphere_volume(const Real &radius);
Real sphere_surface_area(const Real &radius);
void vector_residuals(const std::vector<Real> &vector1, const std::vector<Real> &vector2, Real &max_relative_residual,Real &average_residual);
Real vector_max(const std::vector<Real> &vector);
std::vector<std::string> split(const std::string &str,const std::string &sep);

#endif /* MATERIALENUM_H */

