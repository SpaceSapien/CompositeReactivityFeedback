/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CylindricalMicroCell.cpp
 * Author: chris
 * 
 * Created on April 4, 2017, 10:09 PM
 */

#include "CylindricalMicroCell.h"
#include "CompositeMicroCell.h"
#include <iostream>
#include <cmath>

using namespace MaterialLibrary;

void CylindricalMicroCell::setInnerBoundaryCondition(MicroCellBoundaryCondition* boundary_condition)
{
    delete this->_inner_boundary_condition;
    this->_inner_boundary_condition = boundary_condition;
}

CylindricalMicroCell::CylindricalMicroCell(FuelPinReactor* reactor, const Real &initial_temperature) : MicroCell(reactor)
{
    this->_reactor = reactor;    
    
    //Solver Settings
    _time_step = _reactor->_input_file_reader->getInputFileParameter("Thermal Time Iteration", static_cast<Real>(100e-9) );
    _start_time = 0;
    _outward_integrated_power = 0;
    _integrated_power = 0;
    _current_time = _start_time;    
    
    //Mesh Setup
    _number_mesh_nodes = _reactor->_input_file_reader->getInputFileParameter("Macro Thermal Mesh Size", static_cast<int>(100) );
    
    
    
    _coolant_channel_radius = _reactor->_input_file_reader->getInputFileParameter("Coolant Channel Radius", static_cast<Real>(0.01) );
    Real fuel_pin_pitch = _reactor->_input_file_reader->getInputFileParameter("Hexagonal Pitch", static_cast<Real>(0.04) );
    Real side_length = fuel_pin_pitch / 2.0 ;
    _eqivalent_outer_radius = side_length * sqrt( 3.0 * sqrt(3.0) / ( 2.0 * M_PI) );
    
    
    this->setMesh(new CylindricalMesh(_coolant_channel_radius, _eqivalent_outer_radius, _number_mesh_nodes, _reactor->_number_macro_cells));
    
    _solution.resize( _mesh->numberOfNodes(), initial_temperature );
    
    //Boundary Conditions  
    _outer_boundary_condition = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();   
    _inner_boundary_condition = MicroCellBoundaryCondition::getReflectedHeatFluxBoundaryConditionFactory();   
    
    
    for( int macro_cell_index = 0; macro_cell_index < _reactor->_number_macro_cells; ++macro_cell_index)
    {
        _micro_scale_solvers.push_back(nullptr);
    }
    
}


void CylindricalMicroCell::initializeMicroScaleCells(const std::vector<Real> &macro_cell_power_density)
{
    switch(_reactor->_dimensionality)
    {
        case FuelPinReactor::HomogenousNeutronics :
        case FuelPinReactor::FullHeterogeneous :
        {
            for( int macro_cell_index = 0; macro_cell_index < _reactor->_number_macro_cells; ++macro_cell_index)
            {
                Real temperature, volume;
                
                delete _micro_scale_solvers[macro_cell_index];
                
                this->getAverageCellTemperature( 1, _reactor->_number_macro_cells, (macro_cell_index + 1), temperature, volume);  
                Real macro_scale_radius = _coolant_channel_radius + ( _coolant_channel_radius - _eqivalent_outer_radius ) * ( macro_cell_index + 0.5 );
                
                
                
                Real r_inner = _reactor->_micro_sphere_geometry->_geometry[_number_zones - 2].second;
                Real r_outer = _reactor->_micro_sphere_geometry->_geometry[_number_zones - 1].second * pow(6.0/M_PI, 1.0/3.0);
                Real r_avg_temp = 1.0 / ( 1.0/r_inner - 3.0/( pow(r_outer,3) - pow(r_inner,3)) * ( pow(r_outer,3)/(3*r_inner) - r_outer*r_outer/2.0 + r_inner*r_inner / 6.0 ));
                Real cell_volume = (M_PI*4.0/3.0) * pow(r_outer, 3);
                Real initial_power_density = macro_cell_power_density[macro_cell_index];
                Real kernel_power = cell_volume * initial_power_density;
                Materials matrix_material = _reactor->_micro_sphere_geometry->_geometry[_number_zones - 1].first;
                MaterialDataPacket packet = MaterialLibrary::getMaterialProperties(matrix_material, temperature);
                Real delta_t = kernel_power / (4.0 * M_PI * packet._thermal_conductivity) * ( 1/r_avg_temp - 1/r_outer); 

                        
                FuelPinCompositeMicroCell* new_cell = new FuelPinCompositeMicroCell(_reactor, temperature - delta_t, macro_scale_radius);
                new_cell->setBoundaryCondition(MicroCellBoundaryCondition::getFixedTemperatureBoundaryConditionFactory(temperature - delta_t));
                _micro_scale_solvers[macro_cell_index] = new_cell;                               
                
                std::vector<Real> representative_power_distribution = new_cell->getRespresentativePowerDistribution(initial_power_density);                
                new_cell->iterateInitialConditions(representative_power_distribution);                
            }
            break;
        }
        case FuelPinReactor::HomogenousNeutronicsAndHeatTransfer :
        {
            break;
        }
    }
}


CylindricalMicroCell::~CylindricalMicroCell() 
{
    delete this->_inner_boundary_condition;
    
    //delete the microscale heat transfer cells
    for(auto composite_cell : _micro_scale_solvers)
    {
        delete composite_cell;
    }
}

MicroSolution CylindricalMicroCell::solve(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
{
   switch(_reactor->_dimensionality)
    {
        //Heterogeneous MicroCells
        case FuelPinReactor::HomogenousNeutronics :
        case FuelPinReactor::FullHeterogeneous :
        {
            return solveHeterogeneous(simulation_time_step,power_distribution);
            break;
        }
        //Thermal model is homogenized no need for microcells
        case FuelPinReactor::HomogenousNeutronicsAndHeatTransfer :
        {
            return solveHomogenous(simulation_time_step,power_distribution);
            break;
        }
    }
}

MicroSolution CylindricalMicroCell::solveHeterogeneous(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
{
    std::vector<Real> radial_mesh = _mesh->_position;
    std::vector<Real> previous_solution = this->_solution;
    std::vector<Real> current_solution;
    
    current_solution.resize( previous_solution.size() );
    
   
    
    
    Real time = _current_time;
    unsigned long iteration = 0;
    
    
    while( time < _current_time + simulation_time_step )
    {
       
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        Real dSolution = 0;
        Real outward_heat_flux = 0;
        Real inward_heat_flux = 0; 
        Real net_heat_flux = 0;
        
        //The matrix material is the last material in the geometry
        Materials matrix_material = _reactor->_micro_sphere_geometry->_geometry.back().first;
        std::vector<Real> micro_scale_power = std::vector<Real>(_reactor->_number_macro_cells, -1);
        
        //Important properties of the microcells
        
        for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
        {
            int current_cell = _mesh->_cell_in_zone[radial_index];
            
            //Get some macro scale properties
            Real temperature = previous_solution[radial_index];                        
            MaterialDataPacket material_data = MaterialLibrary::getMaterialProperties(matrix_material, temperature);
            Dimension center_radius = _mesh->getNodeLocation(radial_index);         
            
            Real microcell_internal_power_density = power_distribution[radial_index];
            
            
            //If we haven't calculated out micro to macro cell heat flux calculate it
            if(micro_scale_power[current_cell -1 ] == -1)
            {
                if(this->_micro_scale_solvers[current_cell-1] != nullptr)
                {
                    Real macrocell_volume, macrocell_temperature;
                    this->getAverageCellTemperature(1,_reactor->_number_macro_cells,current_cell,macrocell_temperature,macrocell_volume);
                    micro_scale_power[current_cell-1] = -1 * this->_micro_scale_solvers[current_cell-1]->solveMicroCellCoupling(macrocell_temperature,microcell_internal_power_density, _time_step);
                }
                else
                {
                    micro_scale_power[current_cell-1] = microcell_internal_power_density;
                }
            }
            
            Real internal_power = micro_scale_power[current_cell-1] * _mesh->_volume[radial_index];
            _integrated_power += microcell_internal_power_density * _time_step;
            
            //Check for boundary conditions
            //If we are at the first node
            if( radial_index == 0 )
            {
                
                Real outward_neighbor_cell_temperature = previous_solution[radial_index +  1];
                Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
                Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
                Real outward_dTdr = outward_delta_T / outward_distance; 

                outward_heat_flux =  -_mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity; 
                
                //use the boundary condition to determine the outward heat flux
                inward_heat_flux = this->_inner_boundary_condition->getHeatFlux(this,outward_heat_flux + internal_power);
                _outward_integrated_power += -inward_heat_flux * _time_step;
                _outward_current_power = -inward_heat_flux;
                
            }  
            // If we are at the last node
            else if( radial_index == _mesh->numberOfNodes() - 1)
            {
                //use the boundary condition to determine the outward heat flux
                outward_heat_flux = this->_outer_boundary_condition->getHeatFlux(this,inward_heat_flux + internal_power);
                _outward_integrated_power += -outward_heat_flux * _time_step;
               
            }            
            else
            {   
                //Neighbor dependent Properties
                Real outward_neighbor_cell_temperature = previous_solution[radial_index +  1];
                Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
                Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
                Real outward_dTdr = outward_delta_T / outward_distance; 
                
                
                outward_heat_flux = - _mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity;                               
            }   
            
            Real base_thermal_inertia = material_data._density * material_data._specific_heat;
            Real thermal_inertia = base_thermal_inertia * _mesh->_volume[radial_index];
            
            net_heat_flux = inward_heat_flux + outward_heat_flux + internal_power;  
            dSolution = net_heat_flux / thermal_inertia;
            current_solution[ radial_index ] = previous_solution[ radial_index ] + _time_step * dSolution;
            
            //for the next cell the inward heat flux is just the negative of the outward for the previous cell
            inward_heat_flux = -outward_heat_flux;
            
        }        
        
        previous_solution = current_solution;
        ++iteration;
        time = time + _time_step;
    }
    
    _current_time = time;
    _solution = current_solution;
    return this->getCurrentMicrosolution();
}

MicroSolution CylindricalMicroCell::solveHomogenous(const Real &simulation_time_step,const std::vector<Real> &power_distribution)
{
    
    std::vector<Real> radial_mesh = _mesh->_position;
    std::vector<Real> previous_solution = this->_solution;
    std::vector<Real> current_solution;
    
    current_solution.resize( previous_solution.size() );
    
    Real time = _current_time;
    unsigned long iteration = 0;
    
    while( time < _current_time + simulation_time_step )
    {
       
        if( iteration % 10000 == 0 )
        {
            std::cout << time << "\n"; 
        }
        
        Real dSolution = 0;
        Real outward_heat_flux = 0;
        Real inward_heat_flux = 0; 
        Real net_heat_flux = 0;
        
        for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
        {
            Real temperature = previous_solution[radial_index];                        
            
            MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getHomogenizedMaterialProperties(temperature);
            Dimension center_radius = _mesh->getNodeLocation(radial_index);         
            Real internal_power = power_distribution[radial_index] * _mesh->_volume[radial_index];
            _integrated_power += internal_power * _time_step;
            
            
            //Check for boundary conditions
            //If we are at the first node
            if( radial_index == 0 )
            {
                
                Real outward_neighbor_cell_temperature = previous_solution[radial_index +  1];
                Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
                Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
                Real outward_dTdr = outward_delta_T / outward_distance; 
                
                
                outward_heat_flux =  -_mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity; 
                
                //use the boundary condition to determine the outward heat flux
                inward_heat_flux = this->_inner_boundary_condition->getHeatFlux(this,outward_heat_flux + internal_power);
                _outward_integrated_power += -outward_heat_flux * _time_step;
                _outward_current_power = -outward_heat_flux;
                
            }  
            // If we are at the last node
            else if( radial_index == _mesh->numberOfNodes() - 1)
            {
                //use the boundary condition to determine the outward heat flux
                outward_heat_flux = this->_outer_boundary_condition->getHeatFlux(this,inward_heat_flux + internal_power);
                _outward_integrated_power += -outward_heat_flux * _time_step;
                _outward_current_power = -outward_heat_flux;
            }            
            else
            {   
                //Neighbor dependent Properties
                Real outward_neighbor_cell_temperature = previous_solution[radial_index +  1];
                Real outward_delta_T = temperature - outward_neighbor_cell_temperature;
                Real outward_distance = _mesh->getNodeLocation(radial_index + 1) - center_radius; 
                Real outward_dTdr = outward_delta_T / outward_distance; 
                
                
                outward_heat_flux = - _mesh->_outer_surface[radial_index] * outward_dTdr * material_data._thermal_conductivity;                               
            }   
            
            Real base_thermal_inertia = material_data._density * material_data._specific_heat;
            Real thermal_inertia = base_thermal_inertia * _mesh->_volume[radial_index];
            
            net_heat_flux = inward_heat_flux + outward_heat_flux + internal_power;  
            dSolution = net_heat_flux / thermal_inertia;
            current_solution[ radial_index ] = previous_solution[ radial_index ] + _time_step * dSolution;
            
            //for the next cell the inward heat flux is just the negative of the outward for the previous cell
            inward_heat_flux = -outward_heat_flux;
            
        }        
        
        previous_solution = current_solution;
        ++iteration;
        time = time + _time_step;
    }
    
    _current_time = time;
    _solution = current_solution;
    return this->getCurrentMicrosolution();
}


std::vector<MicroSolution> CylindricalMicroCell::iterateInitialConditions(const std::vector<Real> &power_distribution)
{
    MicroSolution initial_solution = this->presolveSteadyStateAnalytical(power_distribution);
    
    std::vector<MicroSolution> solutions = std::vector<MicroSolution>();
    solutions.push_back(initial_solution);
    
    Real solution_time_step = _time_step * 50000;
    Real max_residual = 1;
    Real avg_residual = 0;
    Real index = 0;
    Real desired_residual =  _reactor->_input_file_reader->getInputFileParameter("Steady State Temperature Solution Max Residual", static_cast<Real>(0.005) );
    
    
    while( max_residual > desired_residual )
    {
       MicroSolution solution = this->solve(solution_time_step,power_distribution); 
       
       if(index > 0)
       {
            //calculate the maximum change in solution from previous to new 
            vector_residuals(solution._solution, solutions.back()._solution, max_residual, avg_residual);
            std::cout<<"Residual: " << max_residual << " Time Step: " << solution_time_step << std::endl;
       }
       solutions.push_back(solution);
       
       index++;
       
    }
    std::vector<Real> macro_cell_power_distribution = this->getMacroCellPowerDensity(power_distribution);
    this->initializeMicroScaleCells(macro_cell_power_distribution);
    
    
    //Redo the iteration this time with the microcells setup
    max_residual = 1;
    avg_residual = 0;
    index = 0;
    desired_residual = _reactor->_input_file_reader->getInputFileParameter("Steady State Temperature Solution Max Residual", static_cast<Real>(0.005) );
    
    
    while( max_residual > desired_residual )
    {
       MicroSolution solution = this->solve(solution_time_step,power_distribution); 
       
       if(index > 1)
       {
            //calculate the maximum change in solution from previous to new 
            vector_residuals(solution._solution, solutions.back()._solution, max_residual, avg_residual);
            std::cout<<"Residual: " << max_residual << " Time Step: " << solution_time_step << std::endl;
       }
       solutions.push_back(solution);
       
       index++;
       
    }
       
    return solutions;
}

std::vector<Real> CylindricalMicroCell::getMacroCellPowerDensity(const std::vector<Real> &node_power_distribution) const
{
    std::vector<Real> macro_cell_power_density;    
    macro_cell_power_density.resize(_reactor->_number_macro_cells);
    
    for( int macro_cell_index = 1; macro_cell_index <= _reactor->_number_macro_cells; ++macro_cell_index )
    {
        Real volume = 0;
        Real volume_weighted_power = 0;
        
        for( std::size_t node_index = 0; node_index <  node_power_distribution.size(); ++node_index)
        {
            if( _mesh->_cell_in_zone[node_index] == macro_cell_index )
            {
                Real current_volume = _mesh->_volume[node_index];
                Real current_power = node_power_distribution[node_index];

                volume_weighted_power += current_volume * current_power;
                volume += current_volume;
            }
        }
        
        Real volume_average_power = volume_weighted_power / volume;
        
        macro_cell_power_density[macro_cell_index-1] =  volume_average_power;
    }
    
    return macro_cell_power_density;
}



MicroSolution CylindricalMicroCell::presolveSteadyStateAnalytical(const std::vector<Real> &power_distribution)
{
    
    
    for( int radial_index = 0; radial_index < _mesh->numberOfNodes(); ++radial_index )
    {
        Real temperature_to_the_left;
        
        //Averaged material properties in the cell assuming right side temperature material properties
        if(radial_index == 0)
        {
            temperature_to_the_left =  _solution[radial_index];
        }
        else
        {
            temperature_to_the_left =  _solution[radial_index - 1];
        }
        
        // Use homogenized data for the full cylinder
        MaterialDataPacket material_data = _reactor->_micro_sphere_geometry->getHomogenizedMaterialProperties(temperature_to_the_left);
        
        Real total_power = 0;
        
        int number_nodes = _mesh->numberOfNodes();
        
        for( int power_index = radial_index; power_index < number_nodes; power_index++ )
        {
             Real volume = _mesh->_volume[power_index];
             total_power += volume * power_distribution[power_index];
        }
       
        Real cell_outer_surface_area = _mesh->_outer_surface[radial_index];
        Real outer_surface_heat_flux = total_power / cell_outer_surface_area;        
        Real power_density = power_distribution[radial_index];    
        Dimension right_radial_position = _mesh->_outer_radius[radial_index];
        Dimension left_radial_position = _mesh->_inner_radius[radial_index];
        
        const Real constant_1 =  1.0 / (2.0 * material_data._thermal_conductivity ) * ( total_power / M_PI  + power_density * right_radial_position * right_radial_position);
        //right_radial_position*right_radial_position * (outer_surface_heat_flux / material_data._thermal_conductivity - power_density * right_radial_position / (3.0 * material_data._thermal_conductivity ) );
        const Real constant_2 = temperature_to_the_left + power_density * left_radial_position * left_radial_position / ( 4.0 * material_data._thermal_conductivity) - constant_1 * log(left_radial_position);

        _solution[radial_index] = -power_density * right_radial_position* right_radial_position / ( 4.0 * material_data._thermal_conductivity) + constant_1 * log(right_radial_position) + constant_2; 
        
    }      
    return MicroSolution(_mesh->_position, _solution,0);
       
}

// Simple homogenous power distribution
std::vector<Real> CylindricalMicroCell::getRespresentativePowerDistribution(const Real &average_power_density)
{
    //Allocate the space for the power distribution
    std::vector<Real> power_distribution = std::vector<Real>();
    Dimension number_points = _mesh->numberOfNodes();
    power_distribution.reserve( number_points);    
    
    for( size_t index = 0; index < number_points; ++index)
    {
        power_distribution.push_back( average_power_density );
    }
    
    return power_distribution;
}


void CylindricalMicroCell::getAverageCellTemperature(const int &zone, const int &zone_divisions, const int &current_division, Real &cell_temperature, Real &cell_volume ) const
{   
    //Since we have only one zone
    Dimension material_inner_radius = _mesh->_inner_radius.front();
    Dimension material_outer_radius = _mesh->_outer_radius.back();
    
    Dimension hex_side = material_outer_radius / sqrt( 3 * sqrt(3) / ( 2.0 * M_PI ) );    
    Dimension hex_close_edge = sqrt(3)/2 * hex_side;
    
    Dimension material_delta_radius = hex_close_edge - material_inner_radius;
    
    
    Dimension inner_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division - 1)/static_cast<Real>(zone_divisions) );
    
    Dimension outer_radius;
            
    if( zone_divisions == current_division )
    {
        outer_radius = material_outer_radius;
    }
    else
    {
        outer_radius = material_inner_radius + material_delta_radius * ( static_cast<Real>(current_division )/static_cast<Real>(zone_divisions) );
    }
    
    this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature, cell_volume);
}

/**
 * Check to see if this mesh at this mesh index is inside the outer and inner radius
 * @param mesh_index
 * @param inner_radius
 * @param outer_radius
 * @return 
 */
inline bool CylindricalMicroCell::inMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
{
    return  _mesh->_inner_radius[mesh_index] < outer_radius && _mesh->_outer_radius[mesh_index] > inner_radius;
}

inline Real CylindricalMicroCell::overlappingVolumeWithMeshCell(const int &mesh_index, const Dimension &inner_radius,const Dimension &outer_radius) const
{
    if( inMeshCell( mesh_index, inner_radius, outer_radius) )
    {
        //if it is completely in the cell
        if( _mesh->_inner_radius[mesh_index] < inner_radius && _mesh->_outer_radius[mesh_index] > outer_radius )
        {
            return _mesh->_volume[mesh_index];
        }
        //The cell over engulfs the cell
        else if( _mesh->_inner_radius[mesh_index] > inner_radius && _mesh->_outer_radius[mesh_index] < outer_radius )
        {
            return _mesh->_volume[mesh_index];
        }
        //There is overlap between regions
        else 
        {
            Real inner_r, outer_r;
            
            //There is left side overlap
            if( _mesh->_inner_radius[mesh_index] < inner_radius )
            {
                inner_r = inner_radius;
            }
            else
            {
                inner_r = _mesh->_inner_radius[mesh_index];
            }
            
            //There is right side overlap
            if(  _mesh->_outer_radius[mesh_index] > outer_radius )
            {
                outer_r = outer_radius;
            }
            else
            {
                outer_r = _mesh->_outer_radius[mesh_index];
            }
            
            return (M_PI * (outer_r * outer_r ) - M_PI * ( inner_r * inner_r));
        }
    }
    else
    {
        return 0;
    }
}

void CylindricalMicroCell::getAverageTemperatureInRadaii(const Dimension &inner_radius, const Dimension &outer_radius, Real &cell_temperature, Real &cell_volume) const
{
    cell_volume = M_PI * outer_radius * outer_radius - M_PI * inner_radius * inner_radius;    
    Real integrated_volume = 0;
    Real temperature_weighted_volume = 0;
    
    for(int index = 0; index < _mesh->numberOfNodes(); index++  )
    {
        Real unit_mesh_volume =  overlappingVolumeWithMeshCell( index, inner_radius, outer_radius);
        temperature_weighted_volume += _solution[index] * unit_mesh_volume;
        integrated_volume += unit_mesh_volume;
    }
    
    cell_temperature = temperature_weighted_volume/cell_volume;
    
}

void CylindricalMicroCell::getAverageZoneTemperature(const int &zone, Real &cell_temperature, Real &cell_volume) const
{
    //For the cylindrical reactor there's only one zone so get the range over the mesh
    Dimension inner_radius = _mesh->_inner_radius.front();
    Dimension outer_radius = _mesh->_outer_radius.back();
    
    return this->getAverageTemperatureInRadaii(inner_radius, outer_radius, cell_temperature,cell_volume);
}

void CylindricalMicroCell::setMesh(CylindricalMesh* mesh)
{
    this->_mesh = mesh;
    this->MicroCell::setMesh(mesh);
}

