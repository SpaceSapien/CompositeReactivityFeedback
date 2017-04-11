/*
 * The MicroCell is a spherical particle
 * 
 */

/* 
 * File:   MicroHeatTransferExplicitSolver.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 5:46 PM
 */
//#define VARIABLE_MATERIAL_PROPERTIES_DIFFERENTIAL
#include <iostream>
#include <vector>
#include <cmath>
#include "MicroSolution.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "MaterialDataPacket.h"
#include "InputDataFunctions.h"
#include "MicroCellBoundaryCondition.h"
#include "MicroCell.h"
#include "InputFileParser.h"
#include "InfiniteCompositeReactor.h"
#include "SphericalMesh.h"



MicroCell::MicroCell(Reactor* reactor) 
{
    this->_reactor = reactor;
    _number_zones = _reactor->_micro_sphere_geometry->_geometry.size();
}

MicroCell::~MicroCell()
{
    delete _mesh;
    delete _outer_boundary_condition;
}

MicroSolution MicroCell::getCurrentMicrosolution() const
{
    MicroSolution solution = MicroSolution(_mesh->_position,_solution,_current_time);
    return solution;
}

void MicroCell::setBoundaryCondition(MicroCellBoundaryCondition* boundary_condition)
{
    if(boundary_condition != nullptr)
    {
        delete _outer_boundary_condition;
    }
    
    this->_outer_boundary_condition = boundary_condition;
}

std::vector<Real> MicroCell::getTallyBasedRepresentativeKernelPowerDistribution(const std::vector<std::vector<Real>> &tally_cell_zone_data, const Real &average_power_density)
{
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    Dimension number_points = _mesh->numberOfNodes();
    
    power_distribution.reserve( number_points);    
    
    Real total_volume = 0;
    Real proportional_power = 0;
    
    for( size_t index = 0; index < number_points; ++index)
    {
        int cell_in_zone = _mesh->_cell_in_zone[index];
        int zone = _mesh->_zone[index];
        
        Real node_power_density = tally_cell_zone_data[zone - 1][cell_in_zone - 1];
        Real node_volume = _mesh->_volume[index];
        Real node_proportional_power = node_power_density * node_volume;
        power_distribution.push_back(tally_cell_zone_data[zone - 1][cell_in_zone - 1]);  
        proportional_power += node_proportional_power;
        total_volume += node_volume;
    }
    
    const Real power_multiplier = average_power_density * total_volume / proportional_power;
    
    for( size_t index = 0; index < number_points; ++index)
    {
         power_distribution[index] *= power_multiplier;
    }    
    
    return power_distribution;
    
}


Real MicroCell::getUnitVolumeIntegratedPower() const
{
    Real boundary_volume = this->getVolume();
    return _outward_integrated_power / boundary_volume;
}
Real MicroCell::getUnitVolumeIntegratedOutwardPower() const
{
    Real boundary_volume = this->getVolume();
    return _integrated_power / boundary_volume; 
}
Real MicroCell::getUnitVolumeOutwardPower() const
{
    Real boundary_volume = this->getVolume();
    return _outward_current_power / boundary_volume;
}
Real MicroCell::getVolume() const
{
    std::vector<Real> unity_volume_weight = std::vector<Real>(_mesh->numberOfNodes(),1);
    Real boundary_volume = _mesh->getVolumeWeightedQuantity(unity_volume_weight);
    return boundary_volume;
}

void MicroCell::setMesh(Mesh* mesh)
{
    this->_mesh = mesh;
}