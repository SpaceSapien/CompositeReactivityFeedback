/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinReactor.h
 * Author: chris
 *
 * Created on March 22, 2017, 10:26 PM
 */

#ifndef FUELPINREACTOR_H
#define FUELPINREACTOR_H

#include "Reactor.h"
#include "CylindricalMicroCell.h"
#include "FuelPinMonteCarlo.h"
#include "CompositeMicroCell.h"

class CylindricalMicroCell;
class FuelPinMonteCarlo;


class FuelPinReactor : public Reactor
{

public:
    
   
    //Create the thermal heat transfer object 
    CylindricalMicroCell* _thermal_solver;
    FuelPinMonteCarlo* _monte_carlo_model;
    std::vector<std::vector<Real>> _microsolver_average_temperatures;
    
    enum DimensionalTreatment
    {
        HomogenousNeutronicsAndHeatTransfer,
        HomogenousNeutronics,
        FullHeterogeneous
    };
    
    static DimensionalTreatment getDimensionalTreatment(const std::string input_text);

    DimensionalTreatment _dimensionality;
    
    //How many macro cells?
    int _number_macro_cells;
    
    
    FuelPinReactor(const std::string &name);
    
    //virtual void timeIterationInnerLoop();
    //virtual void temperatureIterationInnerLoop();
    virtual void initializeReactorProblem();
    virtual void worthStudy();
    
    virtual void setThermalSolver(CylindricalMicroCell* solver);
    virtual void setMoteCarloModel(FuelPinMonteCarlo* mc_model);
    virtual bool significantTemperatureDifferenceInMicroCell(std::vector<Real> comparison);
    virtual bool significantTemperatureDifference(MicroSolution* comparison);
    virtual void monteCarloTimeStepSimulationDataProcessing();
    virtual std::vector<Real> getMicroCellTemperatures();
   
    
    void solveForSteadyStatePowerDistribution(const std::vector<Real> &homogenous_power_density, const Real &initial_power_density);
    void saveMicroScaleDataToFile();
    
    
    virtual ~FuelPinReactor();
    
private:

};

#endif /* FUELPINREACTOR_H */

