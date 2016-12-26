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

SimulationResults::SimulationResults(const std::string &file_name, const std::string &directory,const std::time_t &elapsed_time = -1 ) 
{
    this->readOutputFile(file_name, directory);
    this->_elapsed_time = elapsed_time;
}

SimulationResults::SimulationResults(const SimulationResults &results) 
{
    this->_k_eff = results._k_eff;
    this->_k_eff_sigma = results._k_eff_sigma;
    this->_prompt_neutron_lifetime = results._prompt_neutron_lifetime;
    this->_prompt_neutron_lifetime_sigma = results._prompt_neutron_lifetime_sigma;  
    this->_elapsed_time = results._elapsed_time;
}

SimulationResults::SimulationResults(){}



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
