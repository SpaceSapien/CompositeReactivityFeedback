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
    
    Tally(const Real &value,const Real &uncertainty,const std::string &tally_name);
    Tally(const Real &value,const Real &uncertainty,const std::vector<Real> &energy_bins,const std::vector<Real> &energy_values,const std::vector<Real> &energy_uncertainty,const std::string &tally_name); 
    
    void printTallyInfo(const bool &display_energy_bins = false);
    
    int energyBinSize();
    
    static Tally* getMCNPTally(const std::string &tally_id,const std::string &mcnp_file);
    static std::vector<Real> processMCTALData(const std::string &top_tag,const std::string &tally_string);
    static std::string grabMCNPTallyData(const std::string &tally_id, const std::string &mctal_data, const std::string &mctal_file);
    
private:

};

#endif /* TALLY_H */

