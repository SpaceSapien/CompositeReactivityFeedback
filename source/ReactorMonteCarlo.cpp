/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReactorMonteCarlo.cpp
 * Author: chris
 * 
 * Created on December 17, 2015, 7:09 PM
 */

#include "ReactorMonteCarlo.h"
#include "MicroGeometry.h"
#include <sstream>
#include "InputDataFunctions.h"
#include <chrono>
#include <thread>


ReactorMonteCarlo::ReactorMonteCarlo() {}

ReactorMonteCarlo::ReactorMonteCarlo(MicroCell* &micro_cell_ptr,const Real &starting_k_eff, const  std::string &run_directory)
{
    Real k_eff;
    Real prompt_removal_lifetime;
    Real k_eff_sigma;
    Real prompt_removal_lifetime_sigma;
    
    _run_directory = run_directory;
    
    std::string run_command = "mkdir -p " + _run_directory;
    exec(run_command);
    
    _micro_cell_ptr = micro_cell_ptr;
    getRawCriticalityParameters(k_eff,k_eff_sigma,prompt_removal_lifetime,prompt_removal_lifetime_sigma);
    _virtual_k_eff_multiplier = starting_k_eff / k_eff;  
    _current_k_eff = starting_k_eff;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = k_eff_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
    
    //Parameters that will come from the Monte Carlo Simulation
    std::pair<FissionableIsotope,Real> U235 = { FissionableIsotope::U235, 1 };
    _fission_tally_listing = { U235 };
}

/**
 * Adjusts the k_effective to match the desired starting reactivity
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::updateAdjustedCriticalityParameters()
{
    Real raw_k_effective, raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma;
    
    getRawCriticalityParameters(raw_k_effective,raw_k_effective_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma);    
    
    _current_k_eff = raw_k_effective * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime = prompt_removal_lifetime;
    _current_k_eff_sigma = raw_k_effective_sigma * _virtual_k_eff_multiplier;
    _current_prompt_neutron_lifetime_sigma = prompt_removal_lifetime_sigma;
}

/**
 * Get the raw k_effective from MCNP
 * @param k_eff
 * @param prompt_removal_lifetime
 */
void ReactorMonteCarlo::getRawCriticalityParameters( Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma)
{

    //clean up the previous MCNP data files
    std::string clean_mcnp_command = "rm " + this->_run_directory + "mcnp-composite-fue[l-z].out " + this->_run_directory + "runtp[e-z] " + this->_run_directory + "srct[p-z]";
    exec(clean_mcnp_command);
    
    //create a symbolic link to the Doppler broadened cross sections
    std::string symbolic_link_command = "cd " + this->_run_directory + "; ln -s ../../../doppler-broadened-cs/otf*txt .";
    exec(symbolic_link_command);
    
    //create the MCNP input file
    std::string file_root = "mcnp-composite-fuel";
    std::string input_file_name = file_root + ".inp";
    std::string output_file_name =  file_root + ".out";
    this->createMCNPOutputFile(input_file_name);
   
    std::string command_line_log_file = "mcnp_run_log.txt";
    
    //Run the file
    #ifdef LAPTOP
    
    //We are just running MPI here
    std::string mcnp_path = "/media/chris/DoubleSpace/MCNP/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    std::string command = "cd " + this->_run_directory + "; mpirun -np  7 " + mcnp_path + " i=" + input_file_name + " o=" + output_file_name + " | tee " + command_line_log_file;
    exec(command);
    
    #elif PRACTICE_CLUSTER
    
    std::string mcnp_path = "/share/apps/mcnp/MCNP_CODE/MCNP6/bin/mcnp6.mpi";
    
    //Delete the old log output by overwriting it with a blank string (we are using the log to check to see if the qsub is done)
    std::string remove_command = "cd " + this->_run_directory + "; echo \"\" > " + command_line_log_file;
    exec(remove_command);

    //Run Submission Script with the created MCNP file
    std::string qsub_command = "cd " + this->_run_directory + ";qsub -pe orte 32 ../../../composite-fuel-submission-script.sh " + input_file_name + " " + output_file_name + " " + command_line_log_file;
    exec(qsub_command);
    //Constantly read the output file until it says mcrun done 
    std::string search_lock = "cd " + this->_run_directory + ";cat " + command_line_log_file + " | grep \"mcrun  is done\"";
    std::string is_done;
    
    do
    {   
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        is_done = exec(search_lock);
    }while(is_done == "");
    
    #endif
    

    //Read the output file
    this->readOutputFile(output_file_name, k_eff, k_eff_sigma, prompt_removal_lifetime, prompt_removal_lifetime_sigma);
        
}

void ReactorMonteCarlo::readOutputFile(const std::string &file_name, Real &k_eff, Real &k_eff_sigma, Real &prompt_removal_lifetime, Real &prompt_removal_lifetime_sigma)
{
    
    std::string base_command = "cd " + this->_run_directory + ";cat " + file_name + " | perl -ne ";   
    
    std::string k_eff_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $1'";
    std::string k_eff_command = base_command + k_eff_regex;
    std::string k_eff_str = exec(k_eff_command);
    k_eff = std::stod(k_eff_str);
        
    std::string k_eff_sigma_regex = "'/estimated combined collision\\/absorption\\/track-length keff = ([0-9]\\.[0-9]+) with an estimated standard deviation of ([0-9]\\.[0-9]+)/ && print $2'";
    std::string k_eff_sigma_command = base_command + k_eff_sigma_regex;
    std::string k_eff_sigma_str = exec(k_eff_sigma_command);
    k_eff_sigma = std::stod(k_eff_sigma_str);
    
    std::string prompt_lifetime_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $1'";
    std::string prompt_lifetime_command = base_command + prompt_lifetime_regex;
    std::string prompt_removal_lifetime_str = exec(prompt_lifetime_command);    
    prompt_removal_lifetime = std::stod(prompt_removal_lifetime_str);
        
    std::string prompt_lifetime_sigma_regex = "'/the final combined \\(col\\/abs\\/tl\\) prompt removal lifetime = ([0-9]+\\.[0-9]+E-?[0-9]+) seconds with an estimated standard deviation of ([0-9]+\\.[0-9]+E-?[0-9]+)/ && print $2'";
    std::string prompt_lifetime_sigma_command = base_command + prompt_lifetime_sigma_regex;
    std::string prompt_removal_lifetime_sigma_str = exec(prompt_lifetime_sigma_command);
    prompt_removal_lifetime_sigma = std::stod(prompt_removal_lifetime_sigma_str);
    
    
    /*std::ifstream mcnp_output_file;
    mcnp_output_file.open(file_name);
    
    std::string output_file_text;
    
    if (mcnp_output_file.is_open()) 
    {
        output_file_text.assign( std::istreambuf_iterator<char>(mcnp_output_file) , std::istreambuf_iterator<char>()  );
    }
    mcnp_output_file.close();
    
    
   
    //"estimated combined collision/absorption/track-length keff = 1.79363 with an estimated standard deviation of 0.00064"
    //the final combined (col/abs/tl) prompt removal lifetime = 2.6650E-04 seconds with an estimated standard deviation of 4.6932E-07 */
}

void ReactorMonteCarlo::createMCNPOutputFile(const std::string &file_name)
{
    const Real MeVperK = 8.617e-11;
    Real matrix_temperature = _micro_cell_ptr->getAverageTemperature(1);
    Real fuel_kernel_temperature = _micro_cell_ptr->getAverageTemperature(0);
    std::stringstream mcnp_file;
    Dimension kernel_radius = _micro_cell_ptr->_geometry.getFuelKernelRadius()*100.0; //Dimensions in cm
    Dimension matrix_radius = _micro_cell_ptr->_geometry.getOuterRadius()*100.0;
    
    #ifdef LAPTOP
    
    std::string U238_cs = "92238.66c";
    std::string U235_cs = "92235.66c";
    
    #elif  PRACTICE_CLUSTER 
    
    std::string U238_cs = "92238.80c";
    std::string U235_cs = "92235.80c";    
    
    #endif
     
    mcnp_file << "Composite Fuel Kernel Scale" << std::endl;
    mcnp_file << "c Simulating a small UO2 fuel kernel inside a graphite matrix" << std::endl;
    mcnp_file << "c TMP [MeV] = 8.617e-11 [MeV/K] * T [K]   " << std::endl;
    mcnp_file << "c ----------------------CELL CARDS----------------------------" << std::endl;
    mcnp_file << " 1 1  -10.96  -1      imp:n=1 TMP=" << fuel_kernel_temperature*MeVperK  << " $fuel T = " << fuel_kernel_temperature << " K" << std::endl;
    mcnp_file << " 2 2  -1.78    1 -2   imp:n=1 TMP=" << matrix_temperature*MeVperK <<" $matrix T = " << matrix_temperature << " K" << std::endl;
    mcnp_file << "c end cell cards" << std::endl;
    mcnp_file << std::endl;
    mcnp_file << "c ----------------------SURFACE CARDS---------------------" << std::endl;
    mcnp_file << " 1 SPH 0 0 0 " << kernel_radius << std::endl;
    mcnp_file << "*2 SPH 0 0 0 " << matrix_radius << std::endl;
    mcnp_file << "c end surface cards" << std::endl;
    mcnp_file << std::endl;
    mcnp_file << "c ------------------MATERIAL AND DATA CARDS--------------------" << std::endl;
    mcnp_file << " m1  8016        2          $UO2" << std::endl;
    mcnp_file << "     " << U235_cs << "   0.2" << std::endl;
    mcnp_file << "     " << U238_cs << "   0.8" << std::endl;
    mcnp_file << " mt1 o2-u.27t           $S(a,b) UO2 @ 1200 K" << std::endl;
    mcnp_file << "     u-o2.27t" << std::endl;
    mcnp_file << " m2  6000    1          $Graphite" << std::endl;
    mcnp_file << " mt2 grph.22t           $Graphite S(a,b) treatment @ 500 K" << std::endl;
    mcnp_file << " OTFDB " << U238_cs << std::endl;
    mcnp_file << "       " << U235_cs << std::endl;
    mcnp_file << "       8016.60c" << std::endl;
    mcnp_file << "       6000.60c" << std::endl;
    mcnp_file << " KCODE 5000 1.5 3 33  $need at least 30 active cycles to print results" << std::endl;
    mcnp_file << " KSRC 0 0 0" << std::endl;
    mcnp_file << " print" << std::endl;
    mcnp_file << "c end data" << std::endl;
    
    std::ofstream mcnp_input_file;
    mcnp_input_file.open( this->_run_directory + file_name,std::ios::out);
    mcnp_input_file << mcnp_file.str();
    mcnp_input_file.close();
    
    std::cout<<mcnp_file.str();
}
