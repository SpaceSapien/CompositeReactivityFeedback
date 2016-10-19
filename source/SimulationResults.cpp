/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SimulationResults.cpp
 * Author: chris
 * 
 * Created on March 1, 2016, 1:44 PM
 */
#include "SimulationResults.h"
#include "EnumsAndFunctions.h"
#include "InputDataFunctions.h"

SimulationResults::SimulationResults(const std::string &file_name, const std::string &directory) 
{
    this->readOutputFile(file_name, directory);
}

SimulationResults::SimulationResults(const SimulationResults &results) 
{
    this->_k_eff = results._k_eff;
    this->_k_eff_sigma = results._k_eff_sigma;
    this->_prompt_neutron_lifetime = results._prompt_neutron_lifetime;
    this->_prompt_neutron_lifetime_sigma = results._prompt_neutron_lifetime_sigma;            
}


SimulationResults::SimulationResults() {}

/**
 * Read the MCNP output file and parse the important paramters
 * @param file_name
 * @param k_eff
 * @param k_eff_sigma
 * @param prompt_removal_lifetime
 * @param prompt_removal_lifetime_sigma
 */
void SimulationResults::readOutputFile(const std::string &file_name, const std::string &directory)
{
    
    std::string base_command = "cd " + directory + ";cat " + file_name + " | perl -ne ";   
    
    std::string k_eff_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $1'";
    std::string k_eff_command = base_command + k_eff_regex;
    std::string k_eff_str = exec(k_eff_command);
    _k_eff = std::stod(k_eff_str);
        
    std::string k_eff_sigma_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $2'";
    std::string k_eff_sigma_command = base_command + k_eff_sigma_regex;
    std::string k_eff_sigma_str = exec(k_eff_sigma_command);
    _k_eff_sigma = std::stod(k_eff_sigma_str);
    
    std::string prompt_lifetime_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $1'";
    std::string prompt_lifetime_command = base_command + prompt_lifetime_regex;
    std::string prompt_removal_lifetime_str = exec(prompt_lifetime_command);    
    _prompt_neutron_lifetime = std::stod(prompt_removal_lifetime_str);
        
    std::string prompt_lifetime_sigma_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $2'";
    std::string prompt_lifetime_sigma_command = base_command + prompt_lifetime_sigma_regex;
    std::string prompt_removal_lifetime_sigma_str = exec(prompt_lifetime_sigma_command);
    _prompt_neutron_lifetime_sigma = std::stod(prompt_removal_lifetime_sigma_str);
}

BetaSimulationResults::BetaSimulationResults(const SimulationResults &with_delayed_neutrons, const SimulationResults &without_delayed_neutrons) 
{
    _with_delayed_neutrons = SimulationResults(with_delayed_neutrons );
    _without_delayed_neutrons = SimulationResults(without_delayed_neutrons);
    
    _k_eff = _with_delayed_neutrons._k_eff;
    _k_eff_sigma = _with_delayed_neutrons._k_eff_sigma;
    _prompt_neutron_lifetime = _with_delayed_neutrons._prompt_neutron_lifetime;
    _prompt_neutron_lifetime = _with_delayed_neutrons._prompt_neutron_lifetime_sigma;
    
    _beta_sigma = getBetaEffSigma(_k_eff,_k_eff_sigma, _without_delayed_neutrons._k_eff, _without_delayed_neutrons._k_eff_sigma );
    _beta = ( _k_eff - _without_delayed_neutrons._k_eff ) / _k_eff;
}

/**
 * Helper function to calculate the beta sigma
 * @param k_eff
 * @param k_eff_sigma
 * @param nd_k_eff
 * @param nd_k_eff_sigma
 * @return 
 */
Real BetaSimulationResults::getBetaEffSigma(const Real &k_eff,const Real &k_eff_sigma,const Real &nd_k_eff,const Real &nd_k_eff_sigma)
{
    Real current_beta_eff = ( k_eff - nd_k_eff ) / k_eff;
    
    //sqrt of the sum of the squares of sigmas
    Real difference_uncertainty = sqrt( k_eff_sigma * k_eff_sigma + nd_k_eff_sigma * nd_k_eff_sigma );
    
    //Find the uncertainty fractions
    Real numerator_uncertainty_fraction = difference_uncertainty / ( k_eff - nd_k_eff );
    Real denominator_uncertainty_fraction = k_eff_sigma / k_eff;
    
    //uncertainty is the sum of the squares of the fractions
    return current_beta_eff * sqrt( numerator_uncertainty_fraction * numerator_uncertainty_fraction  + denominator_uncertainty_fraction * denominator_uncertainty_fraction );
}

TallyResults::TallyResults(const std::string &file_name, const std::string &directory) 
{
    this->readTallyFile(file_name, directory);
}




TallyResults::TallyResults() {}

/**
 * Read the MCNP output file and parse the important paramters
 * @param file_name
 * @param k_eff
 * @param k_eff_sigma
 * @param prompt_removal_lifetime
 * @param prompt_removal_lifetime_sigma
 */
void TallyResults::readTallyFile(const std::string &file_name, const std::string &directory)
{
    
   // std::string file_text = this->grabFileVector(directory + file_name);
    
    //std::size_t ntal_start = file_text.find("ntal");
    //std::size_t = file_text.substr(file_text.).find( )
            
            
    /*
    std::string start_index = "tally   " + std::to_string(tally_number);
    start_index.substr()
    std::size_t found = str.find(str2);
    
    if (found!=std::string::npos)
    {
        std::cout << "first 'needle' found at: " << found << '\n';
    }
    
    found=str.find("needles are small",found+1,6);
    
    if (found!=std::string::npos)
    {
        std::cout << "second 'needle' found at: " << found << '\n';
    }
    
    found=str.find("haystack");
    if (found!=std::string::npos)
    {
        std::cout << "'haystack' also found at: " << found << '\n';
    }
    
    found=str.find('.');
    if (found!=std::string::npos)
    {
        std::cout << "Period found at: " << found << '\n';
    }
    // let's replace the first needle:
    str.replace(str.find(str2),str2.length(),"preposition");*/
    
}

