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


class MicroCell
{

public:
    
    MicroCell();
    MicroCell(const MicroGeometry &geometry, const Real &intial_temperature);   
    MicroSolution solve(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution getInitialConditions();
    MaterialDataPacket testMaterialProperties(const Real &radial_position);
    MicroGeometry _geometry;    
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
    Real getAverageTemperature(const int &zone);
    
private:
    
    
    Real calculationMaximumResidual(const std::vector<Real> &vector_1, const std::vector<Real> &vector_2, const Real &time_step);
    MicroSolution solveSecondOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    MicroSolution solveFourthOrder(const Real &simulation_time_step, const std::vector<Real> &power_distribution);
    
    
};

#endif /* MICROHEATTRANSFEREXPLICITSOLVER_H */

