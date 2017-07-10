/* 
 * File:   TallyGroup.cpp
 * Author: chris
 * 
 * Created on November 10, 2016, 9:30 PM
 * Application specific tally file
 */

#include "TallyGroup.h"
#include "Tally.h"
#include <string>
#include "EnumsAndFunctions.h"
#include "InputDataFunctions.h"
#include <vector>
#include <cmath>
#include <fstream>

TallyGroup::TallyGroup(const Real &time,const int &zones,const int &cells_per_zone,
        const std::vector<Tally*> &flux_tallies,
        const std::vector<Tally*> &fission_tallies,
        const std::vector<Tally*> &absorption_rate_tallies,
        const std::vector<Tally*> &fission_rate_tallies,
        const std::vector<Tally*> &capture_rate_tallies
)
{
    _time = time;
    _cells_per_zone = cells_per_zone;
    _flux_tallies = flux_tallies;
    _fission_tallies = fission_tallies;
    _fission_rate_tallies = fission_rate_tallies;
    _absorption_rate_tallies = absorption_rate_tallies;
    _capture_rate_tallies = capture_rate_tallies;
    
    _zones = zones;
}

void TallyGroup::outputTallyToFile(const std::string &file_path, Tally* tally,const std::string &location)
{
    
    std::ofstream output_file;
    output_file.open( file_path, std::ios::app);
    
    output_file << _time << "," << location << "," << tally->_value << "," << tally->_uncertainty;

    for(std::size_t energy_index = 0; energy_index < tally->energyBinSize(); ++energy_index )
    {    
         output_file << "," << tally->_energy_bins[energy_index] << "," << tally->_energy_bin_values[energy_index] << "," << tally->_energy_bin_uncertainties[energy_index];
    }
    
    output_file << std::endl; 
    
     output_file.close();
}

void TallyGroup::outputToFile(const std::string &file_path)
{
    if(file_exists(file_path))
    {
        for( std::size_t index =0; index < this->size(); index++ )
        {
            int current_zone = index / _cells_per_zone + 1;
            int current_cell = index % _cells_per_zone + 1;
            std::string zone_str = std::to_string(current_zone);
            std::string cell_str = std::to_string(current_cell);
            
            
            std::string flux_location =  "Flux-" + zone_str + "-" + cell_str + "," + zone_str + "," + cell_str;
            std::string fission_location =  "Fission-" + zone_str + "-" +cell_str + "," + zone_str + "," + cell_str;
            std::string absorption_rate_location =  "Absorption-Rate-" + zone_str + "-" +cell_str + "," + zone_str + "," + cell_str;
            std::string fission_rate_location =  "Fission-Rate-" + zone_str + "-" +cell_str + "," + zone_str + "," + cell_str;
            std::string capture_rate_location =  "Capture-Rate-" + zone_str + "-" +cell_str + "," + zone_str + "," + cell_str;
            
            
            Tally *flux_tally = _flux_tallies[index];
            Tally *fission_tally = _fission_tallies[index];
            Tally *fission_rate_tally = _fission_rate_tallies[index];
            Tally *absorption_rate_tally = _absorption_rate_tallies[index];
            Tally *capture_rate_tally = _absorption_rate_tallies[index];
            
            this->outputTallyToFile(file_path,flux_tally,flux_location);
            this->outputTallyToFile(file_path,fission_tally, fission_location);            
            this->outputTallyToFile(file_path,absorption_rate_tally, absorption_rate_location);            
            this->outputTallyToFile(file_path,capture_rate_tally, capture_rate_location);            
            this->outputTallyToFile(file_path,fission_rate_tally, fission_rate_location);            
        }    
       
    }
    else
    {
        std::cerr << "Output file " + file_path << " doesn't exist, no tallies are being saved.";
    }
}

int TallyGroup::size()
{
    return _cells_per_zone * _zones;
}

TallyGroup* TallyGroup::MCNPTallyGroupFactory(const std::string &MCTAL_file,const int &zones,const int &cells_per_zone,const Real &time)
{
    int total_tallies = zones * cells_per_zone;
    
    std::vector<Tally*> fission_tallies = std::vector<Tally*>();
    std::vector<Tally*> flux_tallies = std::vector<Tally*>();
    std::vector<Tally*> absorption_rate_tallies = std::vector<Tally*>();
    std::vector<Tally*> fission_rate_tallies = std::vector<Tally*>();
    std::vector<Tally*> capture_rate_tallies = std::vector<Tally*>();
    
    
    for(int index = 1; index <= total_tallies; ++index)
    {
        std::string F7_tally = std::to_string(index) + "7";
        std::string F4_tally = std::to_string(index) + "4";
        
        Tally* fission_tally = Tally::getMCNPTally( F7_tally, MCTAL_file);
        Tally* flux_tally = Tally::getMCNPTally( F4_tally, MCTAL_file);
        
        flux_tallies.push_back(flux_tally);
        fission_tallies.push_back(fission_tally);
    }
    
    //Our Fn4M -1 Cell -6 tally for number of absorptions
    for(int index = total_tallies + 1; index <= total_tallies * 2; ++index)
    {
        std::string F4_absorption_tally = std::to_string(index) + "4";
        
        Tally* absorption_rate_tally = Tally::getMCNPTally( F4_absorption_tally, MCTAL_file);
        
        absorption_rate_tallies.push_back(absorption_rate_tally);
    }
    
    //Our Fn4M -1 Cell -6 tally for number of fissions
    for(int index = total_tallies * 2 + 1; index <= total_tallies * 3; ++index)
    {
        std::string F4_fission_tally = std::to_string(index) + "4";
        
        Tally* fission_rate_tally = Tally::getMCNPTally( F4_fission_tally, MCTAL_file);
        
        fission_rate_tallies.push_back(fission_rate_tally);
    }
    
    //Our Fn4M -1 Cell -6 tally for number of fissions
    for(int index = total_tallies * 3 + 1; index <= total_tallies * 4; ++index)
    {
        std::string F4_capture_tally = std::to_string(index) + "4";
        
        Tally* capture_rate_tally = Tally::getMCNPTally( F4_capture_tally, MCTAL_file);
        
        capture_rate_tallies.push_back(capture_rate_tally);
    }
    
    
    
    return new TallyGroup(time,zones,cells_per_zone, flux_tallies, fission_tallies,absorption_rate_tallies, fission_rate_tallies, capture_rate_tallies);
}

void TallyGroup::print(const bool &print_energies)
{
    for(std::size_t index = 0; index < _flux_tallies.size(); ++index )
    {
        _flux_tallies[index]->printTallyInfo(print_energies);
    }
    for(std::size_t index = 0; index < _fission_tallies.size(); ++index )
    {
        _fission_tallies[index]->printTallyInfo(print_energies);
    } 
    for(std::size_t index = 0; index < _fission_rate_tallies.size(); ++index )
    {
        _fission_rate_tallies[index]->printTallyInfo(print_energies);
    } 
    for(std::size_t index = 0; index < _absorption_rate_tallies.size(); ++index )
    {
        _absorption_rate_tallies[index]->printTallyInfo(print_energies);
    }
    for(std::size_t index = 0; index < _capture_rate_tallies.size(); ++index )
    {
        _capture_rate_tallies[index]->printTallyInfo(print_energies);
    }
}

TallyGroup::~TallyGroup() 
{
    for(std::size_t index=0; index < _flux_tallies.size(); ++index)
    {
        delete _flux_tallies[index];
    }
    
    for(std::size_t index=0; index < _fission_tallies.size(); ++index)
    {
        delete _fission_tallies[index];
    }
    
    for(std::size_t index=0; index < _fission_rate_tallies.size(); ++index)
    {
        delete _fission_rate_tallies[index];
    }
    
    for(std::size_t index=0; index < _absorption_rate_tallies.size(); ++index)
    {
        delete _absorption_rate_tallies[index];
    }
    
    for(std::size_t index=0; index < _capture_rate_tallies.size(); ++index)
    {
        delete _capture_rate_tallies[index];
    }
}

