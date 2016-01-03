/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Solution.cpp
 * Author: chris
 * 
 * Created on November 26, 2015, 1:15 PM
 */
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include "EnumsAndFunctions.h"
#include "MicroSolution.h"
#include "PythonPlot.h"

const std::string MicroSolution::_x_label = "Radial Position [m]";
const std::string MicroSolution::_y_label = "Temperature [K]";

MicroSolution::MicroSolution(const std::vector<Dimension> &dimension, const std::vector<Real> &solution, const Real &time) 
{
    if( solution.size() != dimension.size() )
    {
        throw Errors::SolutionAndGridSizeDifferent;
    }
    
    this->_solution = solution;
    this->_grid = dimension;
    this->_time = time;
}

int MicroSolution::size()
{
    return this->_grid.size();
}

void MicroSolution::plot(const std::string &save_file_name)
{
    std::string title = "Time = " + std::to_string(_time);    
    PythonPlot::plotData(_grid,_solution,_x_label, _y_label,"", title, save_file_name);
}

/**
 * Plots a graph of multiple solutions overlaid on each other
 * 
 * @param plot_data_vector
 * @param number_plots  How many plots to show
 * @param save_file_name
 */
void MicroSolution::plotSolutions(const std::vector<MicroSolution> &plot_data_vector, const int &number_plots, const std::string &save_file_name)
{
    std::vector<std::vector<Real>> dimension_data;
    std::vector<std::vector<Real>> temperature_data;
    std::vector<std::string> legend_data;
    
    
    if(number_plots == 0 || number_plots > plot_data_vector.size())
    {
    
        for( int i_index = 0; i_index < plot_data_vector.size(); i_index++ )
        {
            MicroSolution plot = plot_data_vector[i_index];  

            dimension_data.push_back(plot._grid);
            temperature_data.push_back(plot._solution);
            legend_data.push_back(std::to_string(plot._time));
        }
    }
    else
    {
        int indicies = floor(plot_data_vector.size()/2.0);
        int two_power_index = 1;
        
        while(number_plots > indicies)
        {
            indicies /= 2;
            two_power_index++;
        }
        
        for( int i_index = 0; i_index < plot_data_vector.size(); i_index+=pow(2,two_power_index) )
        {
            MicroSolution plot = plot_data_vector[i_index];  

            dimension_data.push_back(plot._grid);
            temperature_data.push_back(plot._solution);
            legend_data.push_back(std::to_string(plot._time));
        }
    }
    
    
    PythonPlot::plotData(dimension_data,temperature_data, MicroSolution::_x_label, MicroSolution::_y_label,legend_data, "Time History", save_file_name );   
}

void MicroSolution::saveSolutions(const std::vector<MicroSolution> &plot_data_vector, const std::string &save_folder,const std::string &save_file_name)
{
    
    std::ofstream output_file;
    output_file.open( save_folder + save_file_name,std::ios::out);
    
    size_t grid_size = plot_data_vector[0]._grid.size();
    size_t  time_steps= plot_data_vector.size();
    
    //Todo
    //If the file exists append to it else create it
    //if()
    {
        for( size_t i_index = 0; i_index < grid_size; i_index++ )
        {
            output_file << plot_data_vector[0]._grid[i_index] << ", ";
        }
        
        output_file << std::endl;
    }
        
    
    for( size_t i_index = 0; i_index < time_steps; i_index++ )
    {
        MicroSolution plot = plot_data_vector[i_index];  

        for( size_t j_index = 0; j_index < grid_size; j_index++)
        {
            output_file << plot._solution[j_index] << ", ";
        }
        
        output_file << std::endl;
    }
    
    output_file.close();
}