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
#include "MicroGeometry.h"
#include "ExplicitSolverSettings.h"
#include "MicroSolution.h"
#include "MicroCellBoundaryCondition.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"

//forward declare
class InfiniteCompositeReactor;

class MicroCell
{

public:
    
    MicroCell();
    MicroCell(InfiniteCompositeReactor* reactor, const Real &intial_temperature);
    
    InfiniteCompositeReactor* _reactor;
    
    
    MicroSolution getCurrentMicrosolution();
    MicroSolution solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MaterialDataPacket testMaterialProperties(const Real &radial_position);
    ExplicitSolverSettings _solver_settings;
    Real _current_time;
    std::vector<Real> _solution; 
    std::vector<Real> _power_distribution; 
    MicroCellBoundaryCondition _outer_boundary_condition;
    std::vector<MicroSolution> iterateInitialConditions(const Real &initial_average_power_density);
    
    std::vector<Real> getRespresentativePowerDistribution(const Real &average_power_density);
    std::vector<Real> getRespresentativeHomogenizedPowerDistribution(const Real &average_power_density);
    std::vector<Real> getRepresentativeKernelPowerDistribution(const Real &average_power_density);

    void setBoundaryCondition(const MicroCellBoundaryCondition &boundary_condition);
    Real getOuterDerivative();
    void getAverageTemperature(const int &zone, Real &cell_temperature, Real &cell_volume);
    void getCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume );
    
private:
    
    void getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume);
    Real calculationMaximumResidual(const std::vector<Real> &vector_1, const std::vector<Real> &vector_2, const Real &time_step);
    MicroSolution solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    
    
};

#endif /* MICROHEATTRANSFEREXPLICITSOLVER_H */

