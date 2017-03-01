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

WorthStudy::WorthStudy(InfiniteCompositeReactor *reactor) 
{
    _reactor = reactor;
    _output_file = "worth.csv";
}

void WorthStudy::startStudy(const bool &vary_fuel_temperature,const bool &vary_matrix_temperature,const std::string &output_file_name)
{
    //Default geometry sizes
    Real default_sphere_outer_radius = 2e-3;  //meters
    Real default_fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { default_fuel_kernel_outer_radius, default_sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    //Grab geometry from input file
    std::vector<Dimension> dimensions = _reactor->_input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =  _reactor->_input_file_reader->getInputFileParameter("Materials", default_materials);
    
    std::pair<Real,Real> default_worth_temperature_range = { 300, 4000 };
    std::pair<Real,Real> worth_temperature_range= _reactor->_input_file_reader->getInputFileParameter("Worth Range", default_worth_temperature_range);
    Real starting_temperature = worth_temperature_range.first;
    Real ending_temperature = worth_temperature_range.second;
    
    int default_divisions = 12;
    int divisions =  _reactor->_input_file_reader->getInputFileParameter("Worth Temperature Divisions", default_divisions);
    
    Real temperature_increment = (ending_temperature - starting_temperature)/(divisions-1);
    Real current_temperature = starting_temperature;
    
    _micro_sphere_geometry = new MicroGeometry(materials, dimensions);  
    _reactor->_micro_sphere_geometry = _micro_sphere_geometry;
    _monte_carlo_model = new ReactorMonteCarlo(_reactor, _reactor->_results_directory + "run/"); 
    _monte_carlo_model->_tally_cells = true;
    
    for( int division_counter = 0; division_counter < divisions; ++division_counter)
    {
        std::cout << current_temperature << std::endl;
        
        if( vary_fuel_temperature)
        {
            _thermal_solver = new MicroCell(_reactor, current_temperature);
        }
        else
        {
            _thermal_solver = new MicroCell(_reactor, starting_temperature);
        }
        
        if( vary_matrix_temperature )
        {
            _thermal_solver->setOuterMaterialTemperature(current_temperature);
        }
        else
        {
            _thermal_solver->setOuterMaterialTemperature(starting_temperature);
        }
        
        
        
        _reactor->_thermal_solver = _thermal_solver;
        SimulationResults results = _monte_carlo_model->getRawKeff();          
        this->log(results, output_file_name);
        delete _thermal_solver;
        _reactor->_thermal_solver = nullptr;
        
        current_temperature += temperature_increment;
    }
    
    delete _micro_sphere_geometry;
    _reactor->_micro_sphere_geometry = nullptr;
    delete _monte_carlo_model;
}

void WorthStudy::createOutputFile(const std::string &output_file_path)
{
    std::ofstream output_file;
    output_file.open( output_file_path, std::ios::out);    
    output_file << "Fuel Temperature [K],Non Fissile Temperature [K],Power Peaking,K-eigenvalue,K-eigenvalue Sigma,Prompt Neutron Lifetime [s],Prompt Neutron Lifetime Sigma [s],Elapsed Time [s],MC Execution Time [s],Time Per Particle [ms],Time Per Particle CPU [ms/cpu]";
    output_file << std::endl;    
    output_file.close();
}

void WorthStudy::log(const SimulationResults &results,std::string output_file_name)
{
    
    if( output_file_name == "")
    {
        output_file_name = _output_file;
    }
    
    std::ofstream output_file;
    std::string file_path = _reactor->_results_directory + output_file_name;
    
    if(!file_exists(file_path))
    {
        this->createOutputFile(file_path);
    }
    
    output_file.open( file_path, std::ios::app);
    
    //time elapsed since the simulation started
    std::time_t elapsed_time_since_start =  std::time(nullptr) - InfiniteCompositeReactor::_simulation_start_time;   
    Real power_peaking = -1;
    
    if(_monte_carlo_model->_tally_cells)
    {
        std::vector<std::vector<Real>> zones = this->_monte_carlo_model->getZoneCellRelativePowerDensity();
        Real max = zones[0][zones[0].size() - 1];
        Real min = zones[0][0];
        power_peaking = min/max;
    }
    
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

