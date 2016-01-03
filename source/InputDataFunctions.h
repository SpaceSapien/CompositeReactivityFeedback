/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputDataFunctions.h
 * Author: chris
 *
 * Created on November 23, 2015, 10:33 PM
 */


#ifndef INPUTDATAFUNCTIONS_H
#define INPUTDATAFUNCTIONS_H
#include <string>
#include <vector>
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "ReactorKinetics.h"
#include "ReactorMonteCarlo.h"
#include "InfiniteCompositeReactor.h"

std::vector<Real> getInitialTemperatures(std::vector<Dimension> radial_points);
std::vector<Real> getPowerDistribution(std::vector<Dimension> radial_points,const Real kernel_radius,const Real total_power_density);
void plotDelayedPrecursors(const std::vector<Real> &time_record,const std::vector<std::vector<Real>> &delayed_record);
std::string exec(const char* cmd);


#endif /* INPUTDATAFUNCTIONS_H */

