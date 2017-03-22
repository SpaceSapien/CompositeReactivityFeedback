/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MicroGeometry.cpp
 * Author: chris
 * 
 * Created on November 20, 2015, 3:13 PM
 */
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "MicroGeometry.h"
#include "MaterialLibrary.h"
#include "Mixture.h"

using namespace MaterialLibrary;

MicroGeometry::MicroGeometry(const std::vector<Materials> &materials,const std::vector<Dimension> &dimensions) 
{
    
    //Make sure the number of material equals the number of dimensions
    if( ! materials.size() == dimensions.size() )
    {
        throw Errors::IncorrectMaterialsDimensions;
    }
    
    this->_geometry = std::vector< std::pair<Materials,Dimension> >();
    
    
    // Check to make sure that the dimensions are increasing, if so create the 
    // material and dimensions vector
    Dimension last_dimension = 0;
    
    for( int index = 0; index < materials.size(); ++index )
    {
        if( ! dimensions[ index ] > last_dimension )
        {
            throw Errors::NonIncreasingDimensions;
        }
        last_dimension = dimensions[ index ];
        
        std::pair<Materials,Dimension> boundary = { materials[index], dimensions[ index] };
        this->_geometry.push_back( boundary );        
        
    }   
}

MicroGeometry::MicroGeometry()
{
    
}

Dimension MicroGeometry::getOuterRadius() const
{
    return this->_geometry.back().second;
}

Dimension MicroGeometry::getFuelKernelRadius() const
{
    return this->_geometry.front().second;
}

void MicroGeometry::printGeometry()
{
    for( int index = 0; index < this->_geometry.size(); index++ )
    {
        std::cout<< index << "  " 
                 << getMaterialName( _geometry[index].first ) 
                 << "  " << _geometry[index].second << "\n";
    }
}

MaterialDataPacket MicroGeometry::getMaterialProperties(const Real &r,const Real &T)
{
    Materials material = this->getMaterial(r);
      
    return  MaterialLibrary::getMaterialProperties(material,T);
}

Materials MicroGeometry::getMaterial(const Real &r)
{
    Dimension location;
    int number_zones = _geometry.size();
    
    for(int index = 0; index < number_zones; index ++ )
    {
        if( r <= _geometry[index].second )
        {
            return _geometry[index].first;
        }
    }    
    
    throw Errors::rGreaterThanR;
    
}

std::vector<Real> MicroGeometry::getVolumeFractions() const
{
    std::vector<Real> volume_fractions;
    Dimension last_dimension = 0; 
    
    std::size_t number_zones = _geometry.size();
    Dimension outermost_dimension = _geometry[number_zones - 1].second;
    Real total_volume = pow( ( outermost_dimension * 2.0 ) , 3);            
    
    for( std::size_t geometry_index = 0; geometry_index < number_zones; ++geometry_index )
    {  
        
        std::pair<Materials,Dimension> current_material_data = _geometry[geometry_index];
        Dimension current_dimension = current_material_data.second;
        Materials current_material = current_material_data.first;
        Real inner_volume = 0;
        Real outer_volume = 0;
        
        
        inner_volume = sphere_volume(last_dimension);            
        
        //Our outermost shell is cubic
        if( geometry_index == number_zones - 1 )
        {
            //cube volume
            outer_volume = total_volume;
            
        }
        //otherwise it is spherical
        else
        {
            //sphere volume
            outer_volume = sphere_volume(current_dimension);
        }
        
        
        Real current_volume_fraction = (outer_volume - inner_volume)/total_volume;
        
        volume_fractions.push_back(current_volume_fraction);
        last_dimension = current_dimension;
    }
    
    return volume_fractions;
}

MaterialDataPacket MicroGeometry::getHomogenizedMaterialProperties(const Real &T)
{
    
    Real total_density = 0;
    Real total_specific_heat_numerator = 0;    
    Real total_thermal_conductivity = 0;    
    
    Real total_density_dir = 0;
    Real total_specific_heat_numerator_dir = 0;
    Real total_thermal_conductivity_dir = 0;
    
    std::size_t number_zones = _geometry.size();
    
    std::vector<Real> volume_fractions = this->getVolumeFractions();
    
    for( std::size_t geometry_index = 0; geometry_index < number_zones; ++geometry_index )
    {  
        Real current_volume_fraction = volume_fractions[geometry_index];
        
        std::pair<Materials,Dimension> current_material_data = _geometry[geometry_index];
        Materials current_material = current_material_data.first;
        MaterialDataPacket current_material_data_packet = MaterialLibrary::getMaterialProperties(current_material,T);
        
        total_density += current_material_data_packet._density *  current_volume_fraction;
        total_thermal_conductivity += current_material_data_packet._thermal_conductivity * current_volume_fraction;
        total_specific_heat_numerator += current_volume_fraction * current_material_data_packet._density * current_material_data_packet._specific_heat;
        
        total_thermal_conductivity_dir += current_material_data_packet._thermal_conductivity_temperature_derivative * current_volume_fraction;
        total_density_dir += current_material_data_packet._density_temperature_derivative *  current_volume_fraction;
        total_specific_heat_numerator_dir += current_volume_fraction * current_material_data_packet._density * current_material_data_packet._specific_heat_temperature_derivative;
        
    }
    
    Real total_specific_heat = total_specific_heat_numerator / total_density;
    Real total_specific_heat_dir = total_specific_heat_numerator_dir / total_density;
    
    return MaterialDataPacket( total_thermal_conductivity,total_density, total_specific_heat,  total_thermal_conductivity_dir,total_density_dir, total_specific_heat_dir);
}

void MicroGeometry::getHomogenizedMcnpMaterialCard(const int &zone, const std::vector<Real> &geometry_temperature_data, const Real &enrichment_fraction,std::string &material_card_entry,std::string &doppler_card_entry, std::string &mt_card_entry) const
{
    
    std::stringstream material_cards, db_cards, mt_cards;
    
    std::vector<Real> volume_fractions = this->getVolumeFractions();
    
    Mixture mixture;
    
    for( std::size_t volume_index = 0; volume_index < volume_fractions.size(); ++volume_index)
    {
        Compound compound = Compound::getStandardCompound(this->_geometry[volume_index].first, enrichment_fraction);
        Real volume_fraction = volume_fractions[volume_index];        
        mixture.addCompound(compound,volume_fraction);
    }    
    
    Composition<Isotope> pre_modified_isotopic_data = mixture.getIsotopeData();
    std::map<std::string,std::pair<Real,std::string> > isotopic_data = getMCNPCrossSectionLibraries(pre_modified_isotopic_data);
    
    std::vector<std::string> cross_section_library_list;
    
    int isotope_index = 0;
    
    for( auto isotope_data : isotopic_data)
    {
        std::string current_library = isotope_data.first;
        Real mol_fraction = isotope_data.second.first;
        std::string element_name = isotope_data.second.second;
        cross_section_library_list.push_back(current_library);
        
        if(isotope_index == 0)
        {
            material_cards << " m" << std::left << std::setw(3) << zone << "   ";
        }
        else
        {
            material_cards << "        ";
        }
        
        material_cards << std::left << std::setw(14) << current_library << "   " << std::setw(12) << mol_fraction << "     $ " << element_name << std::endl;
        
        ++isotope_index;
    }
    
    std::map<std::string,bool> doppler_broadened_cs = MaterialLibrary::getMcnpDopplerBroadenedCS(cross_section_library_list);
    
    int db_index = 0;
    
    for( auto cross_section : doppler_broadened_cs)
    {
        if(db_index == 0)
        {
            db_cards << " OTFDB ";
        }
        else
        {
            db_cards << "       ";
        }
        
        db_cards << cross_section.first << std::endl;
        
        ++db_index;
    }
    
    mt_cards << " mt" << zone << std::endl;
    
    std::size_t number_materials = _geometry.size();
    
    for(std::size_t material_index = 0; material_index < number_materials; ++material_index)
    {
        Materials current_material = _geometry[material_index].first;
        std::string current_mt_card;
        MaterialLibrary::getMcnpMTCard(current_material, geometry_temperature_data[material_index], current_mt_card);
        mt_cards << current_mt_card;
    }
    
    mt_card_entry = mt_cards.str();
    doppler_card_entry = db_cards.str();
    material_card_entry = material_cards.str();
    
}
