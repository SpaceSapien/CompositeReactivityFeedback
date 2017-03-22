/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   WorthStudy.cpp
 * Author: chris
 * 
 * Created on December 17, 2016, 7:47 PM
 */

#include "WorthStudy.h"
#include "Reactor.h"
#include "CompositeMicroCell.h"
#include "EnumsAndFunctions.h"
#include "InputDataFunctions.h"

WorthStudy::WorthStudy(Reactor *reactor) 
{
    _reactor = reactor;
    _output_file = "worth.csv";
}


void WorthStudy::setThermalSolver(MicroCell* solver)
{
    this->_thermal_solver = solver;
}

void WorthStudy::setMonteCarloModel(ReactorMonteCarlo* monte_carlo_model)
{
    _monte_carlo_model = monte_carlo_model;
}


void WorthStudy::logData(const SimulationResults &results,std::string output_file_name)
{
    
    if( output_file_name == "")
    {
        output_file_name = _output_file;
    }
    
    std::ofstream output_file;
    std::string file_path = _reactor->_results_directory + output_file_name;
    
    if(!file_exists(file_path))
    {
        std::ofstream output_file;
        output_file.open( file_path, std::ios::out);    
        output_file << "Fuel Temperature [K],Non Fissile Temperature [K],Power Peaking,K-eigenvalue,K-eigenvalue Sigma,Prompt Neutron Lifetime [s],Prompt Neutron Lifetime Sigma [s],Elapsed Time [s],MC Execution Time [s],Time Per Particle [ms],Time Per Particle CPU [ms/cpu]";
        output_file << std::endl;    
        output_file.close();
    }
    
    output_file.open( file_path, std::ios::app);
    
    //time elapsed since the simulation started
    std::time_t elapsed_time_since_start =  std::time(nullptr) - InfiniteCompositeReactor::_simulation_start_time;   
    Real power_peaking = this->_monte_carlo_model->getPowerPeaking();
    
    
    
    //Get the fuel and matrix temperature 
    Real fuel_temperature = this->_thermal_solver->_solution.front();
    Real matrix_temperature = this->_thermal_solver->_solution.back();
    std::time_t elapsed_mc_time = results._elapsed_time;
    Real time_per_particle = 1000.0*static_cast<Real>(elapsed_mc_time) / static_cast<Real>( _monte_carlo_model->_k_eff_number_particles); 
    Real time_per_particle_cpu = time_per_particle * this->_monte_carlo_model->_number_cpus;
    
    output_file << fuel_temperature << "," << matrix_temperature << "," << power_peaking << "," << results._k_eff << "," 
                << results._k_eff_sigma << "," << results._prompt_neutron_lifetime << "," 
                << results._prompt_neutron_lifetime_sigma << "," << elapsed_time_since_start << "," << elapsed_mc_time << "," 
                << time_per_particle << "," << time_per_particle_cpu;

    output_file << std::endl;    
    output_file.close();
}