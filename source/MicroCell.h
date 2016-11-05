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
#include "RadialMesh.h"
#include "MicroGeometry.h"
#include "MicroSolution.h"
#include "MicroCellBoundaryCondition.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"

//forward declare for pointers
class MicroCellBoundaryCondition;


class MicroCell
{
    
public:

    enum SolverOrder
    {
        SECOND,
        FOURTH
    };

    
    MicroCell();
    MicroCell(InfiniteCompositeReactor* reactor, const Real &intial_temperature);
    ~MicroCell();
    
    InfiniteCompositeReactor* _reactor;
    RadialMesh* _mesh;
    
    MicroSolution getCurrentMicrosolution();
    MicroSolution solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MaterialDataPacket testMaterialProperties(const Real &radial_position);
    
    Real _time_step;
    Real _start_time;
    int _number_mesh_nodes;
    SolverOrder _solver_order;
    
    
    Real _current_time;
    std::vector<Real> _solution; 
    std::vector<Real> _power_distribution; 
    MicroCellBoundaryCondition* _outer_boundary_condition;
    std::vector<MicroSolution> iterateInitialConditions(const Real &initial_average_power_density);
    
    std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density);
    std::vector<Real> getRespresentativeHomogenizedPowerDistribution(const Real &average_power_density);
    std::vector<Real> getRepresentativeKernelPowerDistribution(const Real &average_power_density);

    void setBoundaryCondition(MicroCellBoundaryCondition* boundary_condition);
    void getAverageTemperature(const int &zone, Real &cell_temperature, Real &cell_volume);
    void getCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume );

    
private:
    
    void getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume);
    Real calculationMaximumResidual(const std::vector<Real> &vector_1, const std::vector<Real> &vector_2, const Real &time_step);
    MicroSolution solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution presolveSteadyStateAnalytical(const Real &average_power_density);
    Real getCellInternalPower(const int &cell, const Real &volumetric_power);
    inline bool inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
    inline Real overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const;
};




#endif /* MICROHEATTRANSFEREXPLICITSOLVER_H */

