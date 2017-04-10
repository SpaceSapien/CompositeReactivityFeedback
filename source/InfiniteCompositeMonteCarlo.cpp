/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InfiniteCompositeMonteCarlo.cpp
 * Author: chris
 * 
 * Created on March 17, 2017, 9:08 AM
 */

#include "InfiniteCompositeMonteCarlo.h"
#include "CompositeMicroCell.h"
#include "MaterialLibrary.h"

using namespace MaterialLibrary;

InfiniteCompositeMonteCarlo::InfiniteCompositeMonteCarlo(InfiniteCompositeReactor* reactor, const std::string &run_dir) : ReactorMonteCarlo(reactor, run_dir) 
{
    this->_reactor = reactor;
}



InfiniteCompositeMonteCarlo::~InfiniteCompositeMonteCarlo() {}

std::string InfiniteCompositeMonteCarlo::getSurfaceCards()
{
    std::stringstream surface_cards;
    
    int surface_card_number = 1;
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    
    //std::string shape_type;        
    
    for( int index = 0; index < _number_zones  ; index++ )
    {
        int current_zone = index + 1;       
        Real last_zone_shell_radius = 0;

        if( index > 0)
        {
            last_zone_shell_radius = geometry_data[index -1 ].second * 100;
        }

        Real zone_shell_radius = geometry_data[index].second * 100; 

        for( int current_surface = 1; current_surface <= _cells_per_zone; current_surface++ )
        {
           Real surface_radius = last_zone_shell_radius + (zone_shell_radius - last_zone_shell_radius ) * (static_cast<Real>(current_surface)/static_cast<Real>(_cells_per_zone));
            
            //we want to put our reflecting boundary condition here
           if( current_zone == _number_zones && current_surface == _cells_per_zone)
           {
               Real box_edge = surface_radius*2;
               Real half_box_edge = surface_radius;
               surface_cards << "*" << surface_card_number << " BOX ";
               surface_cards << -half_box_edge << " " << -half_box_edge << " " << -half_box_edge << "  ";
               surface_cards << box_edge << " 0 0  ";
               surface_cards << "0 " << box_edge << " 0  ";
               surface_cards << "0 0  " << box_edge << std::endl;
           }
           else
           {
               surface_cards << " " << surface_card_number << " SPH 0 0 0 " << surface_radius << std::endl;
           }
           //cm sphere radius

           
           
           surface_card_number++;         
        }
        
    }
    
    return surface_cards.str();
}

std::string InfiniteCompositeMonteCarlo::getMaterialCards()
{
    std::stringstream material_cards, otfdb_card;
    otfdb_card << " OTFDB ";
        
    Real enrichment_fraction = _reactor->_input_file_reader->getInputFileParameter("Uranium Enrichment Fraction", static_cast<Real>(0.2) );

    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
     
       
    for( size_t index = 0; index < _number_zones  ; index++ )
    {
        size_t current_zone = index + 1;
    
        Materials material = geometry_data[index].first; 
        std::string material_card_entry;
        std::string doppler_card_entry;
        
        Real cell_temperature, cell_volume;
        
        int zone = index + 1;
        _reactor->_thermal_solver->getAverageZoneTemperature(zone, cell_temperature, cell_volume);
        
        MaterialLibrary::getMcnpMaterialCard(material,current_zone,cell_temperature,material_card_entry, doppler_card_entry, enrichment_fraction);
        
        material_cards << material_card_entry;
        otfdb_card << doppler_card_entry;        
    }       
    
    
    return material_cards.str() + otfdb_card.str();
}

std::string InfiniteCompositeMonteCarlo::getSingleCellCard(const Materials &material, const int &current_zone, int &cell_number )
{
    std::stringstream cell_card; 
    
     
    //MCNP reads temperature in MeV so we need a conversion factor for our programs K
    const Real MeVperK = 8.617e-11;
    //For now we will assume that the density stays constant regardless of temperature change to perserve conservation of mass
    Real density_derived_temperature = 400;
    //First is the density, second is the density derivative which we don't need. Divide by 1000 to convert from kg/m^3 to g/cm^3
    Real density = MaterialLibrary::getDensityPair(material,density_derived_temperature,0).first/1000;

    for( int current_cell_in_zone = 1; current_cell_in_zone <= _cells_per_zone; current_cell_in_zone++)
    {
        
        //Gather the temperature for the zone in this cell
        Real temperature;
        Real cell_volume;
        
        _reactor->_thermal_solver->getAverageCellTemperature(current_zone - 1, _cells_per_zone, current_cell_in_zone, temperature, cell_volume );
           
        cell_volume *= (100 * 100 * 100);
        
        //we want to put our reflecting boundary condition here
        if( cell_number == 1 )
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << current_zone << " " << std::setw(9) << -density << "  " << std::setw(3) <<                       "  " << std::setw(4)  << -cell_number << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) << getMaterialName(material) << " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
        }
        else
        {
            cell_card << std::left << " " << std::setw(2) << cell_number << " " << std::setw(2) << current_zone << " " << std::setw(9) << -density << " " << std::setw(3) << ( cell_number -1 )  << " " << std::setw(4)  << -cell_number << " imp:n=1 TMP=" << std::setw(11) << temperature*MeVperK << " VOL=" << std::setw(11) << cell_volume << " $ " << std::setw(6) << getMaterialName(material) << " T = " << std::setw(8) << temperature << " K" << " Volume = " << std::setw(11) << cell_volume << " cm^3 " << std::endl;
        }

        cell_number++;
    }
    return cell_card.str();
}


std::string InfiniteCompositeMonteCarlo::getCellCards()
{
    std::stringstream cell_cards;
    
    std::vector<std::pair<Materials, Dimension> > geometry_data = _reactor->_micro_sphere_geometry->_geometry;
    int number_zones = geometry_data.size();   
    int current_cell_number = 1;
    
    //For each zone create a cell
    for( int index = 0; index < number_zones  ; index++ )
    {
        Materials material = geometry_data[index].first;
        int current_zone = index + 1;        
        //Note that the current_cell_number is passed by reference and its value is changed
        cell_cards << this->getSingleCellCard(material,current_zone, current_cell_number);                    
    }
    
    return cell_cards.str();
}
