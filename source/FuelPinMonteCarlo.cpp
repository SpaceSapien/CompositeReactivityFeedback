/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinMonteCarlo.cpp
 * Author: chris
 * 
 * Created on April 5, 2017, 3:24 AM
 */

#include "FuelPinMonteCarlo.h"
#include "CylindricalMicroCell.h"
#include "MaterialLibrary.h"
#include <iomanip>
#include <cmath>
#include <sstream>

using namespace MaterialLibrary;

FuelPinMonteCarlo::FuelPinMonteCarlo(FuelPinReactor* reactor, const std::string &run_dir) : ReactorMonteCarlo(reactor, run_dir) 
{
    this->_reactor = reactor;
    _number_zones = 1;
    
    _coolant_channel_radius = _reactor->_input_file_reader->getInputFileParameter("Coolant Channel Radius", static_cast<Real>(0.01) ) * 100;
    _fuel_pin_pitch = _reactor->_input_file_reader->getInputFileParameter("Hexagonal Pitch", static_cast<Real>(0.04) ) * 100;
    _height = _fuel_pin_pitch / 2;
    _number_macro_cells = _reactor->_number_macro_cells;
    _cells_per_zone = _number_macro_cells;
}



FuelPinMonteCarlo::~FuelPinMonteCarlo() {}

std::string FuelPinMonteCarlo::getSurfaceCards()
{
    std::stringstream surface_cards;
    
   
    
    
    for( int current_surface = 1; current_surface <= _number_macro_cells + 1; current_surface++ )
    {
        
       Real surface_radius =  _coolant_channel_radius +  ( _fuel_pin_pitch/2 * sqrt(3)/2  - _coolant_channel_radius )  * (static_cast<Real>(current_surface-1)/static_cast<Real>(_number_macro_cells));
       
       
       //we want to put our reflecting boundary condition here on the outermost cell
       if( current_surface == _number_macro_cells + 1)
       {
           surface_cards << "*" << current_surface << " RHP ";
           surface_cards << " 0 0 " << -_height/2 << " 0 0 " << _height  << " "  << sqrt(3)/2.0 * _fuel_pin_pitch/2.0 << std::endl;
       }
       //make the center fuel pin
       else 
       {
           surface_cards << " " << current_surface << " RCC ";
           surface_cards << " 0 0 " << -_height << " 0 0 " << _height * 2 << "  " << surface_radius << std::endl;
       }
       //cm sphere radius

       
    }
    
    return surface_cards.str();
}

std::string FuelPinMonteCarlo::getMaterialCards()
{
    MaterialLibrary::MicroGeometry* geometry_data = _reactor->_micro_sphere_geometry;
    
    //Overall average temperature for the mt card
    Real temperature, cell_volume;
    _reactor->_thermal_solver->getAverageZoneTemperature(1, temperature, cell_volume );
    std::vector<Real> geometry_temperature_data = std::vector<Real>(geometry_data->_geometry.size(), temperature);
   
    std::string material_cards, doppler_cards, mt_cards;
    
    Real enrichment_fraction = _reactor->_input_file_reader->getInputFileParameter("Uranium Enrichment Fraction", static_cast<Real>(0.2) );
    //Temperature needed only for mt card
    geometry_data->getHomogenizedMcnpMaterialCard(1, geometry_temperature_data, enrichment_fraction, material_cards, doppler_cards, mt_cards);
      
    
    material_cards += " m2     2004.80c        1\n";

    
    return material_cards + mt_cards + doppler_cards;
    
}

std::string FuelPinMonteCarlo::getCellCards()
{
    std::stringstream cell_card;     
     
    //MCNP reads temperature in MeV so we need a conversion factor for our programs K
    const Real MeVperK = 8.617e-11;
    //For now we will assume that the density stays constant regardless of temperature change to perserve conservation of mass
    Real density_derived_temperature = 400;
    //First is the density, second is the density derivative which we don't need. Divide by 1000 to convert from kg/m^3 to g/cm^3
    
    MaterialDataPacket packet = _reactor->_micro_sphere_geometry->getHomogenizedMaterialProperties(density_derived_temperature);    
    Real density = packet._density/1000;

    for( int cell_number = 1; cell_number <= _number_macro_cells; ++cell_number)
    {
        //Gather the temperature for the zone in this cell
        Real temperature;
        Real cell_volume;
        
        
        if( ! _reactor->_dimensionality == FuelPinReactor::HomogenousNeutronicsAndHeatTransfer )
        {   
            Real microcell_volume, macrocell_temperature;
            _reactor->_thermal_solver->_micro_scale_solvers[cell_number-1]->getAverageZoneTemperature(1, temperature, microcell_volume );
            _reactor->_thermal_solver->getAverageCellTemperature(1, _number_macro_cells, cell_number, macrocell_temperature, cell_volume );
        }
        else 
        {
            _reactor->_thermal_solver->getAverageCellTemperature(1, _number_macro_cells, cell_number, temperature, cell_volume );
        }
        
        
           
        cell_volume *= (100.0 * 100.0 * 100.0) * (  _height / 100 ); // divide by 100 because its in meters
        
        if( cell_number != _number_macro_cells )
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << 1 << " " << std::setw(9) << -density << " " << std::setw(3) << ( cell_number )  << " " << std::setw(4)  << -(cell_number + 1) << -(_number_macro_cells + 1) << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) <<  " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
        }
        else
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << 1 << " " << std::setw(9) << -density << " " << std::setw(3) << ( cell_number )  << " " << std::setw(4)  << -(_number_macro_cells + 1) << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) <<  " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
        }
    }
    
    Real coolant_channel_volume = M_PI * _coolant_channel_radius * _coolant_channel_radius * _height;
    
    cell_card << std::left << " " << (_number_macro_cells + 1) << " " << std::setw(2) << 2 << " -0.0178" << "  -1 " << -(_number_macro_cells + 1) << " imp:n=1 TMP=" << 800*MeVperK  << " VOL=" << coolant_channel_volume << std::endl;
    
    return cell_card.str();
}

/*std::vector< std::vector<Real> > FuelPinMonteCarlo::getZoneCellRelativePowerDensity()
{
    std::vector<std::vector<Real>> cell_zone_power_densities = std::vector<std::vector<Real>>();
    
    Real maximum_power_density = -1;
    
    //Go through each zone and each cell in each zone 
    for(int zone = 1; zone <= _number_zones; zone++)
    {
        std::vector<Real> cell_power_densities = std::vector<Real>();

        for( int cell = 1; cell <= _number_macro_cells; ++cell)
        {   
            Real power_density;
            
            //If no tallies set up so that the first zone is powered
            if(! _tally_cells || _tally_groups.size() == 0)
            {
                
                if( cell == 1 && zone == 1)
                {
                    power_density = 1.0;
                }
                else
                {
                    power_density = 0.0;
                }
            }
            //Otherwise use tallies
            else
            {
                TallyGroup* latest_tally_group = this->_tally_groups.back();    
                power_density = latest_tally_group->_fission_tallies[ (zone-1)*_number_macro_cells + (cell-1) ]->_value;               
            }
            
            if(power_density > maximum_power_density)
            {
                maximum_power_density = power_density;
            }
            
            cell_power_densities.push_back(power_density);            
        }
        
        cell_zone_power_densities.push_back(cell_power_densities);
    }
    
    for(int zone = 0; zone < _number_zones; ++zone)
    {
        for( int cell = 0; cell < _number_macro_cells; ++cell)
        {               
            cell_zone_power_densities[zone][cell] = cell_zone_power_densities[zone][cell]/maximum_power_density;
        }
     }
    
    return cell_zone_power_densities;
}*/