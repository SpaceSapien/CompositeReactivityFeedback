/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMicroCell.h
 * Author: chris
 *
 * Created on March 15, 2017, 7:05 AM
 */

#ifndef HOMOGENOUSMICROCELL_H
#define HOMOGENOUSMICROCELL_H
#include "MicroCell.h"
#include "InfiniteHomogenousReactor.h"
#include "HomogenousMesh.h"

class InfiniteHomogenousReactor;

class HomogenousMicroCell : public MicroCell
{
public:

    HomogenousMesh* _mesh;
    
    HomogenousMicroCell(InfiniteHomogenousReactor* reactor, const Real &initial_temperature);
    virtual ~HomogenousMicroCell();    
    virtual MicroSolution solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    virtual std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density);
    
    virtual void setMesh(HomogenousMesh* mesh);
    virtual void getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const;
    virtual void getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const;
    virtual Real getVolume() const;

    
private:
    
    

};

#endif /* HOMOGENOUSMICROCELL_H */

