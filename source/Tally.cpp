/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Tally.cpp
 * Author: chris
 * 
 * Created on March 5, 2016, 12:51 AM
 */

#include <string>
#include <regex>
#include <iostream>
#include <iomanip>
#include "Tally.h"
#include "EnumsAndFunctions.h"
#include "InputDataFunctions.h"



Tally::Tally(const Real &value, const Real &uncertainty,const std::string &tally_name) : _value(value), _uncertainty(uncertainty), _tally_name(tally_name)
{
    _energy_bins = std::vector<Real>();
    _energy_bin_values = std::vector<Real>();
    _energy_bin_uncertainties = std::vector<Real>();
}

Tally::Tally(const Real &value,const Real &uncertainty,const std::vector<Real> &energy_bins,const std::vector<Real> &energy_values,const std::vector<Real> &energy_uncertainty,const std::string &tally_name) :
_value(value), _uncertainty(uncertainty), _energy_bins(energy_bins), _energy_bin_values(energy_values), _energy_bin_uncertainties(energy_uncertainty), _tally_name(tally_name)
{
    
}

std::string Tally::grabMCNPTallyData(const std::string &tally_id, const std::string &mctal_data,const std::string &mcnp_mctal_file)
{
    std::vector<std::string> split_data = split(mctal_data, "\n");
    
    
    bool found_tally = false;
    std::size_t start_index = 0;
    std::size_t end_index = 0;
    
    for( std::size_t index=0; index < split_data.size(); ++index)
    {
        std::string tally_identifier = "tally " + tally_id;
        
        if(split_data[index].find(tally_identifier) != std::string::npos)
        {
            std::cout<<split_data[index];
            start_index = index;
            found_tally = true;
            break;
        }
    }
    
    if(! found_tally)
    {
        throw "Tally F" + tally_id + " doesn't exist in " + mcnp_mctal_file;
    }
    
    bool found_tfc = false;
    
    for( std::size_t index=start_index; index < split_data.size(); ++index)
    {
        if(split_data[index].find("tfc") != std::string::npos)
        {
            end_index = index + 1;
            found_tfc = true;
            break;
        }
    }
    
    if(! found_tfc)
    {
        throw "Tally F" + tally_id + " doesn't have a tfc " + mcnp_mctal_file;
    }
    
    std::string single_tally_data = "";
    
    for( std::size_t index=start_index; index <= end_index; ++index)
    {
        single_tally_data += split_data[index] + "\n";
    }
    
    return single_tally_data;
    
}

Tally* Tally::getMCNPTally(const std::string &tally_id, const std::string &mcnp_mctal_file)
{
    if(! file_exists( mcnp_mctal_file ))
    {
         throw "Tally file doesn't exist: line " +std::to_string( __LINE__ ) + "file" + __FILE__;        
    }
    
    
    std::string file_text = get_file_text(mcnp_mctal_file);
    
    //std::cout<< file_text;
    
    std::regex multiple_spaces("[ ]{2,}");
    std::string compressed_text = std::regex_replace(file_text, multiple_spaces, " ");
    
    //std::cout << compressed_text;
    
    std::smatch match;
    /*std::regex tally_info("tally " + tally_id + "(.|\\r|\\n)+?tfc"); // [-+.0-9E ]+\\n( [-+.0-9E ]+\\n){1,}");
    
    
    if(! std::regex_search(compressed_text, match, tally_info))
    {
        throw "Tally F" + tally_id + " doesn't exist in " + mcnp_mctal_file;
    }*/
    
    std::string tally_mcnp_text = Tally::grabMCNPTallyData(tally_id, compressed_text, mcnp_mctal_file);//match.str();       
    std::regex value_and_uncertainty("tfc.+?\\n [-+.0-9E]+ ([-+.0-9E]+) ([-+.0-9E]+)");
    
    if(! std::regex_search(tally_mcnp_text, match, value_and_uncertainty))
    {
        throw "Tally F" + tally_id + " has no tfc line in " + mcnp_mctal_file;
    }
   
    Real value = static_cast<Real>(std::stod(match.str(1)));
    Real uncertainty = static_cast<Real>(std::stod(match.str(2)));
    
    std::vector<Real> energy_bins = Tally::processMCTALData("et", tally_mcnp_text);
    
    Tally* tally;
    
    if(energy_bins.size() != 0)
    {
        std::vector<Real> tally_values = std::vector<Real>();
        std::vector<Real> tally_uncertainty = std::vector<Real>();
        std::vector<Real> output_values = Tally::processMCTALData("vals", tally_mcnp_text);
                
        for(int index=0; index < energy_bins.size(); index++)
        {
            tally_values.push_back( output_values[ index * 2 ]);
            tally_uncertainty.push_back(output_values[ index * 2 + 1]);
        }
        
        tally = new Tally(value, uncertainty, energy_bins, tally_values, tally_uncertainty, tally_id);
    }
    else
    {
        tally = new Tally(value, uncertainty,tally_id);
    }
    return tally;
    
}

std::vector<Real> Tally::processMCTALData(const std::string &top_tag,const std::string &tally_mcnp_text)
{
    std::regex energy_tally_check( top_tag + "(.+)?\\n");
    std::vector<Real> output = std::vector<Real>();
    std::smatch match;
    
    //If there is an energy tally
    if( std::regex_search( tally_mcnp_text, match, energy_tally_check ) )
    {
        int start_bin_cutoff = match.position() + match.length();
       
        std::regex energy_groups_regex(top_tag + "(.+)?\\n([0-9E+-\\. ]+\\n){1,}");
        std::regex_search(tally_mcnp_text, match, energy_groups_regex);
        int length = match.length();
        int position = match.position();
        int end_bin_cutoff = length + position;
        
        std::string raw_energy_bins = tally_mcnp_text.substr(start_bin_cutoff,end_bin_cutoff-start_bin_cutoff);
        
        std::regex remove_return("\\n");
        std::string energy_bins = std::regex_replace(raw_energy_bins, remove_return, "");
        
        
        std::regex number("([0-9E+-.]+)");
        std::sregex_iterator start_bin = std::sregex_iterator(energy_bins.begin(), energy_bins.end(), number);
        std::sregex_iterator end_bin = std::sregex_iterator();
        
        
        for (std::sregex_iterator i = start_bin; i != end_bin; ++i) 
        {
            std::smatch match = *i;
            std::string match_str = match.str();       
            Real data = static_cast<Real>(std::stod(match_str));
            output.push_back(data);            
        }

        
    }
    return output;
}

void Tally::printTallyInfo(const bool &energy_bins)
{
    std::cout<< std::left << std::setw(8) << _tally_name << " Value: "<< std::setw(11) << _value << "  sigma="<< _uncertainty;
    
   
    std::cout << std::endl;
    
    if( _energy_bins.size() > 0 && energy_bins)
    {
        Real last_bin = 0;
        
        for( std::size_t index = 0; index < _energy_bins.size(); ++index)
        {
            std::cout << "Energy " << std::left << std::setw(12) << last_bin << " to " << std::setw(12) << _energy_bins[index] << " MeV    Values: " << std::setw(12) << _energy_bin_values[index] << " sigma=" << _energy_bin_uncertainties[index] << std::endl;
            last_bin = _energy_bins[index];
        }
    }
    
    
}

int Tally::energyBinSize()
{
    return this->_energy_bin_values.size();
}