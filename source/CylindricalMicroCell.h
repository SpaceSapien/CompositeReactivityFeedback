/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CylindricalMicroCell.h
 * Author: chris
 *
 * Created on April 4, 2017, 10:09 PM
 */

#ifndef CYLINDRICALMICROCELL_H
#define CYLINDRICALMICROCELL_H

#include "FuelPinReactor.h"
#include "MicroSolution.h"
#include "Reactor.h"
#include "MicroCell.h"
#include "CylindricalMesh.h"
#include "FuelPinCompositeMicroCell.h"

class FuelPinReactor;
class FuelPinCompositeMicroCell;

using namespace MaterialLibrary;


class CylindricalMicroCell : public MicroCell
{
public:
    
    std::vector<FuelPinCompositeMicroCell*> _micro_scale_solvers;
    
    
    CylindricalMesh* _mesh;
    MicroCellBoundaryCondition* _inner_boundary_condition;
    FuelPinReactor* _reactor;
    Real _coolant_channel_radius;
    Real _eqivalent_outer_radius;
    
    
    CylindricalMicroCell(FuelPinReactor* reactor, const Real &initial_temperature);
    virtual ~CylindricalMicroCell();
    virtual MicroSolution solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution);
    virtual std::vector<MicroSolution> iterateInitialConditions(const std::vector<Real> &power_distribution);
   
    virtual void setMesh(CylindricalMesh* mesh);
    //If tallies aren't being taken this is the representative power distribution
    std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density);
    
    virtual void setInnerBoundaryCondition(MicroCellBoundaryCondition* boundary_condition);
    
    virtual void getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const;
    virtual void getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const;
 
    void initializeMicroScaleCells(const std::vector<Real> &macro_cell_power_density);
    
protected:
    
    MicroSolution solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    virtual MicroSolution presolveSteadyStateAnalytical(const std::vector<Real> &radial_power_density);
    virtual void getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume) const;
    virtual inline bool inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
    virtual inline Real overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
    std::vector<Real> getMacroCellPowerDensity(const std::vector<Real> &node_power_distribution) const;
    MicroSolution solveHeterogeneous(const Real &simulation_time_step,const std::vector<Real> &power_distribution);
    MicroSolution solveHomogenous(const Real &simulation_time_step,const std::vector<Real> &power_distribution);
    
    

};

#endif /* CYLINDRICALMICROCELL_H */

