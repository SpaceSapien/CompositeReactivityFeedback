/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CompositeMicroCell.h
 * Author: chris
 *
 * Created on March 13, 2017, 1:25 PM
 */

#ifndef COMPOSITEMICROCELL_H
#define COMPOSITEMICROCELL_H
#include "InfiniteCompositeReactor.h"
#include "MicroSolution.h"
#include "Reactor.h"
#include "MicroCell.h"
#include "SphericalMesh.h"


using namespace MaterialLibrary;

class CompositeMicroCell : public MicroCell
{
public:
    
    enum SolverOrder
    {
        SECOND,
        FOURTH
    };
    
    SolverOrder _solver_order;
    SphericalMesh* _mesh;
    
    CompositeMicroCell(InfiniteCompositeReactor* reactor, const Real &initial_temperature);
    virtual ~CompositeMicroCell();
    virtual MicroSolution solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution);
    virtual std::vector<MicroSolution> iterateInitialConditions(const std::vector<Real> &power_distribution);
   
    virtual void setMesh(SphericalMesh* mesh);
    //If tallies aren't being taken this is the representative power distribution
    std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density);
    
    virtual void setOuterMaterialTemperature(const Real &outer_temperature);
    
    virtual void getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const;
    virtual void getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const;

    
protected:
    
    MicroSolution solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    virtual MicroSolution presolveSteadyStateAnalytical(const std::vector<Real> &radial_power_density);
    virtual void getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume) const;
    virtual inline bool inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
    virtual inline Real overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
    
    

};

#endif /* COMPOSITEMICROCELL_H */

