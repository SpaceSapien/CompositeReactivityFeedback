/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Tally.h
 * Author: chris
 *
 * Created on March 5, 2016, 12:51 AM
 */

#ifndef TALLY_H
#define TALLY_H
#include "EnumsAndFunctions.h"

class Tally 
{
public:
    
    std::string _tally_name;
   
    std::vector<Real> _energy_bins;
    std::vector<Real> _energy_bin_uncertainties;
    std::vector<Real> _energy_bin_values;
    Real _value;
    Real _uncertainty;
    
    Tally(Real value, Real uncertainty, std::string tally_name);
    Tally(Real value, Real uncertainty, std::vector<Real> energy_bins, std::vector<Real> energy_values, std::vector<Real> energy_uncertainty, std::string tally_name); 
    
    void printTallyInfo(bool display_energy_bins = false);
    
    int energyBinSize();
    
    static Tally* getMCNPTally(std::string tally_id, std::string mcnp_file);
    static std::vector<Real> processMCTALData(std::string top_tag, std::string tally_string);
    
private:

};

#endif /* TALLY_H */

