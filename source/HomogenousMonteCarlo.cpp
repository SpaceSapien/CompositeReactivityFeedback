/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HomogenousMonteCarlo.cpp
 * Author: chris
 * 
 * Created on March 17, 2017, 12:55 PM
 */

#include "HomogenousMonteCarlo.h"
#include "MaterialDataPacket.h"
#include <sstream>
#include <iomanip>

HomogenousMonteCarlo::HomogenousMonteCarlo(InfiniteHomogenousReactor* reactor, const std::string &run_directory) : ReactorMonteCarlo(reactor, run_directory)
{
    _reactor = reactor;
    //Overrride to make one zone and one cell
    this->_number_zones = 1;
    this->_cells_per_zone = 1;
}



HomogenousMonteCarlo::~HomogenousMonteCarlo() {}

std::string HomogenousMonteCarlo::getMaterialCards()
{
    MaterialLibrary::MicroGeometry* geometry_data = _reactor->_micro_sphere_geometry;
    
    Real temperature, cell_volume;
    _reactor->_thermal_solver->getAverageCellTemperature(1, 1, 1, temperature, cell_volume );
    std::vector<Real> geometry_temperature_data = std::vector<Real>(geometry_data->_geometry.size(), temperature);
   
    std::string material_cards, doppler_cards, mt_cards;
    
    Real enrichment_fraction = _reactor->_input_file_reader->getInputFileParameter("Uranium Enrichment Fraction", static_cast<Real>(0.2) );

    
    geometry_data->getHomogenizedMcnpMaterialCard(1, geometry_temperature_data, enrichment_fraction, material_cards, doppler_cards, mt_cards);
    
    return material_cards + mt_cards + doppler_cards;
}

std::string HomogenousMonteCarlo::getCellCards()
{
    std::stringstream cell_card;
    
    MaterialLibrary::MicroGeometry* geometry_data = _reactor->_micro_sphere_geometry;
    //MCNP reads temperature in MeV so we need a conversion factor for our programs K
    const Real MeVperK = 8.617e-11;
    //For now we will assume that the density stays constant regardless of temperature change to perserve conservation of mass
    Real density_derived_temperature = 400;
    //First is the density, second is the density derivative which we don't need. Divide by 1000 to convert from kg/m^3 to g/cm^3
    MaterialDataPacket packet = geometry_data->getHomogenizedMaterialProperties(density_derived_temperature);
    
    Real temperature,cell_volume;
    _reactor->_thermal_solver->getAverageCellTemperature(1, 1, 1, temperature, cell_volume );
    
    Real density = packet._density/1000; //density in g/cm^3
    cell_volume = 1 * (100 * 100 * 100);  // one cubic meter in cm^3        
    cell_card << std::left << " " << std::setw(2) << 1 << " " << std::setw(2) << 1 << " " << std::setw(9) << -density << "  " << std::setw(3) <<                       "  " << std::setw(4)  << -1 << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) << "Homogenized" << " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
      
    
    return cell_card.str();
}

std::string HomogenousMonteCarlo::getSurfaceCards()
{

    std::stringstream surface_cards;   
    int surface_card_number = 1;    
    // 1 m^3 box with a reflected boundary condition
    surface_cards << "*1 BOX  -50 -50 -50 100 0 0  0 100 0 0 0  100" << std::endl; 
    return surface_cards.str();
}
