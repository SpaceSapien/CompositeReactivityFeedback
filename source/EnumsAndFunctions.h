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

#endif /* MATERIALENUM_H */

