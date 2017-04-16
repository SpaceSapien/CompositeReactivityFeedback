/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Reactor.cpp
 * Author: chris
 * 
 * Created on March 12, 2017, 10:34 PM
 */

#include "Reactor.h"
#include "PythonPlot.h"
#include "EnumsAndFunctions.h"
#include "InputDataFunctions.h"
#include "ReactivityInsertion.h"
#include <string>
#include <iostream>



/** The computer start time for the calculation */
const std::time_t Reactor::_simulation_start_time = std::time(nullptr);


/**
 * Helper function to turn the string of the input file into the enum
 * @param string_type the input file's string for the recalculation type
 * @return enum of the recalculation type
 */
Reactor::MoneCarloRecalculation Reactor::getRecalculationType(const std::string &string_type)
{
    if(string_type == "Temperature")
    {
        return Reactor::MoneCarloRecalculation::Temperature;
    }
    else if( string_type == "Time" )
    {
        return Reactor::MoneCarloRecalculation::Time;
    }
}


/** **Incomplete**
 * This constructor loads the conditions of a previously run program and continues it until the new end time
 * @param old_results_folder  the folder with the old files in it
 * @param new_end_time        the new end time (must be longer than the old one
 */
/*
Reactor::Reactor(const std::string old_results_folder, Real new_end_time) 
{    
    //Check to make sure the old results folder exists 
    if( ! file_exists(old_results_folder) )
    {
        std::cerr << "Results Directory: " + old_results_folder + " Doesn't exist";
        throw 1;

    }

    //Then check for the files that are needed
    std::string old_input_file = old_results_folder + "/input_file.inp";
    std::string old_data_file = old_results_folder + "/datafile.csv";
    std::string old_temeprature_file = old_results_folder + "/temperature-data.csv";

    if( ! file_exists(old_input_file) || ! file_exists(old_data_file) || ! file_exists(old_temeprature_file))
    {
        std::cerr << "Missing Data Files " + old_input_file + " " + old_data_file + old_temeprature_file;
        throw 1;
    }
    
    //Initialize the input file reader
    this->_input_file_reader = new InputFileParser( old_input_file );
    
    //Create the Results Folder
    time_t run_identification_number = std::time(nullptr);
    
    std::string default_run_name = "Unnamed-Run";
    _run_name = _input_file_reader->getInputFileParameter("Run Name", default_run_name);
    
    _results_directory =  "results/" + _run_name + "-" + std::to_string(run_identification_number) + "/";    
    
    _data_file = "datafile.csv";
    std::string folder_command = "mkdir -p " + _results_directory;
    exec( folder_command );
    
    //Copy the files to the new directory
    std::string copy_input_file_command = "cp " + old_input_file + " " + _results_directory + "input_file.inp";
    exec( copy_input_file_command );
    
    std::string copy_old_data_file = "cp " + old_data_file + " " + _results_directory + "input_file.inp";
    exec( copy_old_data_file );
    
    std::string copy_temperature_file_command = "cp " + old_temeprature_file + " " + _results_directory + "input_file.inp";
    exec( copy_temperature_file_command );
    
    //_monte_carlo_number_iterations = -1;
    //_transient_time = -1;
    
    
}
*/

bool Reactor::_otf_sab = false;

/**
 * Setup the directory structure of the problem then run the pre transient calculations
 * @param input_file_name the name and path of the input file
 */
Reactor::Reactor(const std::string &input_file_name) 
{
    //Check to make sure the old dire
    if( ! file_exists(input_file_name) )
    {
        std::cerr << "Input File: " + input_file_name + " doesn't exist";
        throw 1;

    }
    
    //Initialize the input file reader
    this->setInputFileParser(new InputFileParser( input_file_name ));
    
    //Create the Results Folder
    time_t run_identification_number = std::time(nullptr);
    
    std::string default_run_name = "Unnamed-Run";
    _run_name = _input_file_reader->getInputFileParameter("Run Name", default_run_name);    
    _results_directory =  "results/" + _run_name + "-" + std::to_string(run_identification_number) + "/";        
    _data_file = "datafile.csv";    
    
    std::string folder_command = "mkdir -p " + _results_directory;
    exec( folder_command );
    
    //Copy the input file to the run folder
    std::string copy_input_file_command = "cp " + input_file_name + " " + _results_directory + "input_file.inp";
    exec( copy_input_file_command );
    
    Reactor::_otf_sab = _input_file_reader->getInputFileParameter("OTFSAB", false);    
    
    //These will be created at the child class level
    _micro_sphere_geometry= nullptr;
    _thermal_solver = nullptr;
    _kinetics_model = nullptr;
    _monte_carlo_model = nullptr;
    _reactivity_insertion_model = nullptr;
}


Reactor::~Reactor() 
{
    delete _input_file_reader;        
    delete _micro_sphere_geometry;
    delete _thermal_solver;
    delete _kinetics_model;
    delete _monte_carlo_model;
    delete _reactivity_insertion_model;    
}


/**
 * Perform all steps and logic associated with a time based MC eignvalue timestep
 */
void Reactor::timeIterationInnerLoop()
{
    Real last_reported_time = 0;
        
    //This loop iterates the kinetics and thermal model data
    for( _inner_time_step = 0 ; _inner_time_step < _monte_carlo_time_iteration ; _inner_time_step += _kinetics_thermal_sync_time_step)
    {
        //Solve the kinetics model
        Real current_power = _kinetics_model->solveForPower(_kinetics_thermal_sync_time_step);

        //Get a vector representation of the radial power distribution
        std::vector<Real> power_distribition = _thermal_solver->getRespresentativePowerDistribution( current_power /* current shape */ );

        //Get the thermal solution
        _thermal_solver->solve( _kinetics_thermal_sync_time_step, power_distribition); 

        //Add the current time steps solution to the record (we can record all of them every 20 us let usually do every millisecond)            
        if( ( _inner_time_step - last_reported_time )  >= _power_and_delayed_neutron_record_time_step )
        {
            last_reported_time = _inner_time_step;

            Real absolute_time = _transient_time + _inner_time_step;

            //Everytime the MC is recalculated store this data
            std::pair<Real,Real> power_entry = { absolute_time, current_power };
            _power_record.push_back(power_entry);

            std::pair<Real,std::vector<Real>> delayed_entry = { absolute_time, _kinetics_model->_delayed_precursors };
            _delayed_record.push_back(delayed_entry);
        }

    }

    this->monteCarloTimeStepSimulationDataProcessing();
    
    //if there is still enough time left to do another monte carlo time iteration
    if( _transient_time + _inner_time_step < _end_time)
    {
        Real last_k_eff = _reactivity_insertion_model->getCurrentKeff(_transient_time);
        
        _monte_carlo_model->updateAdjustedCriticalityParameters();

        Real current_k_eff =  _reactivity_insertion_model->getCurrentKeff(_transient_time);

        Real difference = current_k_eff - last_k_eff;

        Real k_eff_change = std::abs( difference);


        Real k_eff_sigma = _reactivity_insertion_model->getCurrentKeffSigma(_transient_time);
        Real min_k_eff_change = k_eff_sigma*2/3;
        Real max_k_eff_change = k_eff_sigma*2.5;

        //We need smaller time steps
        if( max_k_eff_change < k_eff_change )
        {
            _monte_carlo_time_iteration *= 0.5;
        }
        //we can get away with bigger time steps
        else if( min_k_eff_change > k_eff_change )
        {
            _monte_carlo_time_iteration *= 1.5;
        }       
    }
    
     _monte_carlo_number_iterations++;
    
}

/**
 * Perform all steps and logic associated with a temperature based MC eignvalue timestep
 */
void Reactor::temperatureIterationInnerLoop()
{
    Real last_reported_time = 0;
    
    std::vector<std::vector<Real>> tally_power_density_map = _monte_carlo_model->getZoneCellRelativePowerDensity();
        
    //This loop iterates the kinetics and thermal model data
    for( _inner_time_step = 0 ; true /*infinite loop this and exit with a break*/ ; _inner_time_step += _kinetics_thermal_sync_time_step)
    {
        //Solve the kinetics model
        Real current_power = _kinetics_model->solveForPower(_kinetics_thermal_sync_time_step);
        std::vector<Real> power_distribition;
        
        Real absolute_time = _transient_time + _inner_time_step;
        
        if(_monte_carlo_model->_tally_cells)
        {
            //Get a vector representation of the radial power distribution
            power_distribition = _thermal_solver->getTallyBasedRepresentativeKernelPowerDistribution(tally_power_density_map, current_power );
        }
        else
        {
            power_distribition = _thermal_solver->getRespresentativePowerDistribution( current_power );
        }
        //Get the thermal solution
        MicroSolution current_solution = _thermal_solver->solve( _kinetics_thermal_sync_time_step, power_distribition); 

        
        //Record the power and delayed records on an intermediate time scale
        if( ( _inner_time_step - last_reported_time )  >= _power_and_delayed_neutron_record_time_step )
        {
            last_reported_time = _inner_time_step;
            
            //Everytime the MC is recalculated store this data
            std::pair<Real,Real> power_entry = { absolute_time, current_power };
            _power_record.push_back(power_entry);

            std::pair<Real,std::vector<Real>> delayed_entry = { absolute_time, _kinetics_model->_delayed_precursors };
            _delayed_record.push_back(delayed_entry);
        }
        
        //If there has been a significant temperature change recalculate the k_eff
        if
        (  
            this->significantTemperatureDifference(&current_solution) || 
            _reactivity_insertion_model->rampNeedsReactivityMonteCarloUpdate(_transient_time , _inner_time_step) ||
            absolute_time > _end_time
        )
        {
            break;
        }
        
        //If this is a ramp insertion and there has been a certain amount of time passed record in the log 
        if(this->_reactivity_insertion_model->rampNeedsReactivityLogUpdate(_transient_time , _inner_time_step) )
        {
            this->monteCarloTimeStepSimulationDataProcessing();
        }

    }

    
    _monte_carlo_model->updateAdjustedCriticalityParameters();
    _monte_carlo_time_iteration = _inner_time_step;   
    _monte_carlo_number_iterations++;
}


/**
 * Start a simulation at transient time zero and continue at inner time steps according to the keignevalue recalculation methods
 */
void Reactor::simulateTransient()
{
    bool worth_study = _input_file_reader->getInputFileParameter("Worth", false);
    
    if(worth_study)
    {
        this->worthStudy();
    }
    
    _monte_carlo_number_iterations = 0;
    _transient_time = 0;
    
    this->initializeReactorProblem();
    
    //Reset the integrated power out and integrated power metrics in preparation for the transient
    this->_thermal_solver->_outward_integrated_power = 0;
    this->_thermal_solver->_integrated_power = 0;
    
    //Simulate the transient the outer loop is the monte carlo simulation
    for( _transient_time = 0; _transient_time < _end_time; _transient_time += _inner_time_step)
    {
        this->monteCarloTimeStepSimulationDataProcessing();
        
        if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Time)
        {
            this->timeIterationInnerLoop();        
        }
        else if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Temperature )
        {
            this->temperatureIterationInnerLoop();
        }
        else
        {
            throw "Monte Carlo Iteration Type Unknown";
        }
    }
    
    //Save data and create the graphs for the post simulation data
    this->postSimulationProcessing();
}


void Reactor::initializeReactorProblem()
{
    _monte_carlo_number_iterations = 0;
    _transient_time = -1;
    
    //Are we recalculating the k-eignvalue on a timestep or on a temperature bases
    //Temperature is the newer better way to do it
    std::string mc_recalc_type = _input_file_reader->getInputFileParameter(std::string("Monte Carlo Recalculation Type"), std::string("Temperature"));
    _monte_carlo_reclaculation_type = InfiniteCompositeReactor::getRecalculationType(mc_recalc_type);
    
    if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Temperature )
    {
        _maximum_allowed_temperature_difference = _input_file_reader->getInputFileParameter("MC Recalc Delta T", static_cast<Real>(20.0));
        _maximum_allowed_average_temperature_difference = _input_file_reader->getInputFileParameter("MC Recalc Avg Delta T", static_cast<Real>(_maximum_allowed_temperature_difference/3.3));
        _monte_carlo_time_iteration = -1;
    }
    else if( _monte_carlo_reclaculation_type == MoneCarloRecalculation::Time )
    {
        _monte_carlo_time_iteration =  _input_file_reader->getInputFileParameter("Monte Carlo Recalculation Timestep", static_cast<Real>(0.01) );  //How often to calculate keff and the prompt neutron lifetime
    }
    
    //This command establishes the logging of python plot commands. Useful for debugging
    PythonPlot::_log_file = this->_results_directory + "graph_log.log";
    
    //Default geometry sizes
    Real default_sphere_outer_radius = 2e-3;  //meters
    Real default_fuel_kernel_outer_radius = 4e-4;
    std::vector<Real> default_dimensions =  { default_fuel_kernel_outer_radius, default_sphere_outer_radius };
    std::vector<Materials> default_materials = { Materials::UO2, Materials::C }; 
    
    //Grab geometry from input file
    std::vector<Dimension> dimensions = _input_file_reader->getInputFileParameter("Radaii", default_dimensions);
    std::vector<Materials> materials =   _input_file_reader->getInputFileParameter("Materials", default_materials);
    
    //Create a new geometry
    this->setMicroSphereGeometry(new MaterialLibrary::MicroGeometry(materials, dimensions));   
    
   
    
    
    //Time stepping parameters
    _power_and_delayed_neutron_record_time_step =  _input_file_reader->getInputFileParameter("Power Record", static_cast<Real>(0.0005) );  //How often to calculate keff and the prompt neutron lifetime
    _kinetics_thermal_sync_time_step = _input_file_reader->getInputFileParameter("Kinetics Thermal Data Sync", static_cast<Real>(20e-6) );      //How often to couple the kinetics and heat transfer routines    
    _end_time = _input_file_reader->getInputFileParameter("Calculation End Time", static_cast<Real>(1.00) );                                    //How many seconds should the simulation last 
  
}

/**
 * Comparing the current microsolution to the one on record as the last k-eignvalue time step microsolution
 * 
 * @param comparison the current microsolution
 * @return true or false if there is a significant enought temperature change to warrant a new k-eignevalue calc
 */
bool Reactor::significantTemperatureDifference(MicroSolution* comparison)
{
     MicroSolution* reference = &_plot_solutions.back();
     
     MicroSolution differences = MicroSolution::temperatureDifference( comparison, reference);
     
     Real max_difference = 0;
     Real average_difference = 0;
     
     for(std::size_t ii=0; ii < differences.size(); ++ii)
     {
         if(max_difference < differences._solution[ii])
         {
             max_difference = differences._solution[ii];
         }
         
         if(ii !=0)
         {
            average_difference += differences._solution[ii] * ( std::pow(differences._grid[ii],3) - std::pow(differences._grid[ii - 1],3) );         
         }
         else
         {
             average_difference += differences._solution[ii] * ( std::pow(differences._grid[ii],3) );         
         }
     }     
     
     average_difference /= std::pow(differences._grid.back(),3);
     
     if( max_difference > _maximum_allowed_temperature_difference || average_difference > _maximum_allowed_average_temperature_difference )
     {
         return true;
     }
     else
     {
         return false;
     }     
}

void Reactor::createOutputFile()
{
    
    std::ofstream output_file;
    output_file.open( this->_results_directory + this->_data_file, std::ios::out);
    
    output_file << "Iteration,Time [s],Timestep [s],Power [W/m^3],k_eff,k_eff sigma,neutron lifetime [s],Neutron Lifetime sigma [s],"
                << "Beta_eff,Beta_eff sigma,Run Time [s],Edge Temp [K],Gamma,Power Peaking,Current Power Out [W/m^3],"
                << "Integrated Outward Power [W*s/m^3],Integrated Power [W*s/m^3],MC Execution Time [s],Time Per Particle [ms],"
                << "Time Per Particle CPU [ms/cpu]";
    
    for(size_t index = 1; index <= 6; index++ )
    {
        output_file << ",Group " << index;
    }
    
    output_file << std::endl;    
    output_file.close();
}

/**
 * Create an output file if not created and then write the current eigenvalue's timestep
 * 
 * @param time
 * @param power
 * @param k_eff
 * @param k_eff_sigma
 * @param neutron_lifetime
 * @param neutron_lifetime_sigma
 * @param beta_eff
 * @param beta_eff_sigma
 * @param hot_temperature
 */
void Reactor::saveCurrentData(const Real &time, const Real &power, const Real &k_eff, const Real &k_eff_sigma, const Real &neutron_lifetime, const Real &neutron_lifetime_sigma, const Real &beta_eff, const Real &beta_eff_sigma, const Real &hot_temperature, const Real &gamma, const Real &integrated_power_out, const Real &integrated_power, const Real &current_power_out)
{
    std::ofstream output_file;
    std::string file_path = this->_results_directory + this->_data_file;
    
    if(!file_exists(file_path))
    {
        this->createOutputFile();
    }
    
    output_file.open( file_path, std::ios::app);
    
    //time elapsed since the simulation started
    std::time_t elapsed_time_since_start =  std::time(nullptr) - InfiniteCompositeReactor::_simulation_start_time;   
    std::time_t mc_execution_time = _monte_carlo_model->_current_mc_exection_elapsed_time;
    
    Real time_per_particle = 1000.0*static_cast<Real>(mc_execution_time) / static_cast<Real>( _monte_carlo_model->_current_number_particles); 
    Real time_per_particle_cpu = time_per_particle * this->_monte_carlo_model->_number_cpus;
    
    Real power_peaking = -1;
    
    if(_monte_carlo_model->_tally_cells)
    {
        std::vector<std::vector<Real>> zones = this->_monte_carlo_model->getZoneCellRelativePowerDensity();
        Real max = zones[0][zones[0].size() - 1];
        Real min = zones[0][0];
        power_peaking = min/max;
    }
    
    output_file << _monte_carlo_number_iterations << "," << time << "," << _monte_carlo_time_iteration << "," << power << "," <<  k_eff << "," << k_eff_sigma
                << "," << neutron_lifetime << "," << neutron_lifetime_sigma << "," << beta_eff << "," << beta_eff_sigma << "," << elapsed_time_since_start 
                << "," << hot_temperature << "," << gamma << "," << power_peaking << "," << current_power_out<< "," << integrated_power_out << "," << integrated_power 
                << "," << mc_execution_time << "," << time_per_particle << "," << time_per_particle_cpu;
    
    auto delayed_precursors = _kinetics_model->_delayed_precursors;
    
    for(size_t index = 0; index < delayed_precursors.size(); index++ )
    {
        output_file << "," << delayed_precursors[index];
    }

    output_file << std::endl;    
    output_file.close();
}


void Reactor::monteCarloTimeStepSimulationDataProcessing()
{
    //Set the current power
    Real current_power = this->_kinetics_model->_current_power;

    //Gather the parameters from the monte carlo model 
    //The Monte Carlo model is run on the outer loop
    Real prompt_removal_lifetime = _monte_carlo_model->_current_prompt_neutron_lifetime;
    Real k_eff = _reactivity_insertion_model->getCurrentKeff(_transient_time); 
    Real k_eff_sigma = _reactivity_insertion_model->getCurrentKeffSigma(_transient_time);
    Real beta_eff = _monte_carlo_model->_current_beta_eff;
    Real beta_eff_sigma = _monte_carlo_model->_current_beta_eff_sigma;        
    Real lambda = prompt_removal_lifetime/k_eff;
    Real lambda_sigma = _monte_carlo_model->_current_prompt_neutron_lifetime_sigma/k_eff;

    //Create the data structures for each time step's storage
    std::tuple<Real,Real,Real> k_effective_data = std::make_tuple( _transient_time , k_eff, k_eff_sigma );
    _k_eff_record.push_back( k_effective_data );

    std::tuple<Real,Real,Real> beta_effective_data = std::make_tuple( _transient_time , beta_eff, beta_eff_sigma );
    _beta_eff_record.push_back( beta_effective_data );

    std::tuple<Real,Real,Real> prompt_removal_lifetime_pair = std::make_tuple(_transient_time, lambda, lambda_sigma);
    _prompt_life_time_record.push_back(prompt_removal_lifetime_pair);

    Real reactivity = (k_eff - 1.0)/k_eff;
    Real reactivity_uncertainty = std::abs( reactivity ) * ( std::abs(k_eff_sigma / ( k_eff - 1)) + std::abs( k_eff_sigma / k_eff ) );
    Real reactivity_pcm = 10000.0 * reactivity;
    Real reactivity_pcm_uncertainty = 10000.0 * reactivity_uncertainty; 
    Real reactivity_cents = reactivity / beta_eff;
    Real reactivity_cents_uncertainty = std::abs( reactivity_cents ) * ( std::abs( beta_eff_sigma / beta_eff) + std::abs( reactivity_uncertainty / reactivity ));

    std::tuple<Real,Real,Real> reactivity_pcm_pair = std::make_tuple(_transient_time, reactivity_pcm, reactivity_pcm_uncertainty);
    _reactivity_pcm_record.push_back(reactivity_pcm_pair);        

    std::tuple<Real,Real,Real> reactivity_cents_pair = std::make_tuple(_transient_time, reactivity_cents, reactivity_cents_uncertainty);
   _reactivity_cents_record.push_back(reactivity_cents_pair);

    Real hottest_temperature = this->_thermal_solver->_solution[0]; 
    Real gamma = this->_reactivity_insertion_model->getCurrentVirtualKeffMultiplier(_transient_time);
    
    
    Real outward_energy_flux = _thermal_solver->getUnitVolumeIntegratedOutwardPower(); 
    Real integrated_power = _thermal_solver->getUnitVolumeIntegratedPower(); 
    Real current_power_out = _thermal_solver ->getUnitVolumeOutwardPower(); 
    
    
    this->saveCurrentData(_transient_time, current_power, k_eff, k_eff_sigma, lambda, lambda_sigma, beta_eff, beta_eff_sigma, hottest_temperature, gamma, outward_energy_flux, integrated_power, current_power_out);
    
    MicroSolution solution = this->_thermal_solver->getCurrentMicrosolution();
    std::vector<MicroSolution> current_solution = { solution };
    _plot_solutions.push_back(solution); 
    MicroSolution::saveSolutions( current_solution, this->_results_directory );
}

void Reactor::postSimulationProcessing()
{
    //Fill in the data logs for the last time step
    this->monteCarloTimeStepSimulationDataProcessing();
    
    //Post Processing graph creation
    MicroSolution::plotSolutions( _plot_solutions, 6 , this->_results_directory + "solutions-graph.png");
    PythonPlot::plotData(      _power_record,            "Time [s]", "Power Density [W/m^3]",      "", "Power vs. Time",                   this->_results_directory + "power-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _prompt_life_time_record, "Time [s]", "Prompt Neutron Lifetime [s]","", "Prompt Neutron Lifetime vs. Time", this->_results_directory + "prompt-neutron-lifetime-graph.png", {0, _end_time} );
    PythonErrorPlot::plotData( _k_eff_record,            "Time [s]", "K effective",    "",             "K-eff vs. Time",                   this->_results_directory + "k-eff-graph.png",                   {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_cents_record, "Time [s]", "Reactivity [Dollars]",    "",    "Reactivity vs. Time",              this->_results_directory + "reactivity-cents-graph.png",        {0, _end_time} );
    PythonErrorPlot::plotData( _reactivity_pcm_record,   "Time [s]", "Reactivity [pcm]",    "",        "Reactivity vs. Time",              this->_results_directory + "reactivity-pcm-graph.png",          {0, _end_time} );
    PythonErrorPlot::plotData( _beta_eff_record,         "Time [s]", "Beta effective",    "",          "Beta-eff vs. Time",                this->_results_directory + "beta-eff-graph.png",                {0, _end_time} );
    PythonPlot::plotData(      _delayed_record,          "Time [s]", "Delayed Precursors",         {}, "Keff vs. Delayed Precursors",      this->_results_directory + "delayed-precursors.png",            {0, _end_time} );
    PythonPlot::createPlots();        
}

void Reactor::setThermalSolver(MicroCell* solver)
{
    this->_thermal_solver = solver;
}

void Reactor::setMicroSphereGeometry(MaterialLibrary::MicroGeometry* geometry)
{
    this->_micro_sphere_geometry = geometry;
}

void Reactor::setKineticsModel(ReactorKinetics* kinetics_model)
{
    this->_kinetics_model = kinetics_model;
}

void Reactor::setMoteCarloModel(ReactorMonteCarlo* monte_carlo)
{
    this->_monte_carlo_model = monte_carlo;
}

void Reactor::setReactivityInsertionModel(ReactivityInsertion* insertion_model)
{
    this->_reactivity_insertion_model = insertion_model;
}

void Reactor::setInputFileParser(InputFileParser* file_reader)
{
    this->_input_file_reader = file_reader;
}