/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MicroHeatTransferExplicitSolver.h
 * Author: chris
 *
 * Created on November 20, 2015, 5:46 PM
 */


#ifndef MICROCELL_H
#define MICROCELL_H
#include <vector>
#include "SphericalMesh.h"
#include "MicroCellBoundaryCondition.h"
#include "MicroSolution.h"
/*#include "MicroGeometry.h"
#include "MicroSolution.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"
*/
//forward declare for pointers
class MicroCellBoundaryCondition;
//class InfiniteCompositeReactor;
class Reactor;


class MicroCell
{
    
public:

    MicroCell(Reactor* reactor);
    ~MicroCell();
    
    Reactor* _reactor;
    Mesh* _mesh;
    
    Real _time_step;
    Real _start_time;
    int _number_mesh_nodes;
    Real _outward_integrated_power;
    Real _outward_current_power;
    Real _integrated_power;

    virtual MicroSolution getCurrentMicrosolution() const;
    virtual MicroSolution solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution) = 0;
 
    virtual Real getUnitVolumeIntegratedPower() const;
    virtual Real getUnitVolumeIntegratedOutwardPower() const;
    virtual Real getUnitVolumeOutwardPower() const;
    virtual Real getVolume() const;
    
    Real _current_time;
    std::vector<Real> _solution; 
    std::vector<Real> _power_distribution; 
    MicroCellBoundaryCondition* _outer_boundary_condition;
    
    virtual void setBoundaryCondition(MicroCellBoundaryCondition* boundary_condition);
    
    virtual void getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const = 0;
    virtual void getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const = 0;

    virtual std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density) = 0;
    std::vector<Real> getTallyBasedRepresentativeKernelPowerDistribution(const std::vector<std::vector<Real>> &tally_cell_zone_data, const Real &average_power_density);
    
    virtual void setMesh(Mesh* mesh);
    
    
};




#endif /* MICROHEATTRANSFEREXPLICITSOLVER_H */

