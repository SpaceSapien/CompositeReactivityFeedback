/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PythonPlot.cpp
 * Author: chris
 * 
 * Created on December 27, 2015, 7:14 PM
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "PythonPlot.h"
#include "InputDataFunctions.h"

PythonPlot::PythonPlot()
{
    
}

PythonPlot::PythonPlot(const std::string &x_data, const std::string &y_data, const std::string &x_label, const std::string &y_label,  const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits )
{
    _legend = legend_data;
    _save_file = save_file_name;
    _title = title_data;
    _x_data = x_data;
    _x_label = x_label;
    _y_data = y_data;
    _y_label = y_label;
    _y_limits = y_limits;
    _x_limits = x_limits;    
}


void PythonPlot::plotData(const std::string& x_data, const std::string& y_data, const std::string& x_label, const std::string& y_label, const std::string& legend_data, const std::string& title_data, const std::string& save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits)
{
    PythonPlot plot = PythonPlot(x_data,y_data,x_label,y_label,legend_data,title_data,save_file_name, x_limits, y_limits);
    plot.plot();
}


void PythonPlot::plotData(const std::vector<std::vector<Real>> &x_data, const std::vector<std::vector<Real>> &y_data, const std::string &x_label, const std::string &y_label,const std::vector<std::string> &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits )
{
    if( x_data.size() == y_data.size() && y_data.size() == legend_data.size() )
    {
        std::string legend_string = "";
        
        for(size_t index = 0; index < x_data.size(); index++ )
        {
            if(x_data[index].size() != y_data[index].size())
            {
                throw Error::XYDataVectorsNotSameSize;
            }
            
            legend_string += legend_data[index] + " ";
        }
        
        std::string x_data_string = PythonPlot::commandLinePlotData(x_data);
        std::string y_data_string = PythonPlot::commandLinePlotData(y_data);
        PythonPlot plot = PythonPlot(x_data_string,y_data_string,x_label,y_label,legend_string,title_data,save_file_name, x_limits, y_limits);
        plot.plot();
    }
    else
    {
        throw Error::XYDataVectorsNotSameSize;
    }
}


void PythonPlot::plotData(const std::vector<std::pair<Real,Real> > &data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits )
{
    std::vector<Real> x_data;
    std::vector<Real> y_data;
    
    for(int x=0; x < data.size(); x++ )
    {
        x_data.push_back(data[x].first);
        y_data.push_back(data[x].second);
    }
    PythonPlot::plotData(x_data,y_data,x_label,y_label,legend_data,title_data,save_file_name, x_limits, y_limits);
       
}

void PythonPlot::plotData(const std::vector<std::pair<Real, std::vector<Real> > > &data, const std::string &x_label, const std::string &y_label,const std::vector<std::string> &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits )
{
 
    std::vector< std::vector<Real> > y_data;
    std::vector< std::vector<Real> >x_data;
    std::vector<Real> temp_x_data;
    std::vector<std::string> legend_data_passing;
    
    //If there is no legend data passed in
    if(legend_data.size() == 0)
    {
        //create the legend and the vector<vector<Real>> of times
        for( size_t group_index = 0; group_index < data[0].second.size(); group_index++)
        {
            std::string legend_string = "Group-" + std::to_string(group_index + 1);
            legend_data_passing.push_back(legend_string);
        }
    }
    else
    {
        legend_data_passing = legend_data;
    }
    
    
    //create the legend and the vector<vector<Real>> of times
    for( size_t group_index = 0; group_index < data[0].second.size(); group_index++)
    {
        std::vector<Real> delayed_group = std::vector<Real>();
        y_data.push_back(delayed_group);        
    }
    
    //We have to flip the vectors
    for( size_t time_index = 0; time_index < data.size(); time_index++ )
    {
        for( size_t group_index = 0; group_index < data[time_index].second.size(); group_index++)
        {
           y_data[group_index].push_back(data[time_index].second[group_index]);
        }
        
        temp_x_data.push_back(data[time_index].first);
    }
    
    
    for(int x=0; x < y_data.size(); x++ )
    {
        x_data.push_back(temp_x_data);
    }
        
    PythonPlot::plotData(x_data, y_data,x_label,y_label,legend_data_passing,title_data,save_file_name, x_limits, y_limits);
}

void PythonPlot::plotData(const std::vector<Real> &x_data, const std::vector<Real> &y_data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits)
{
    if( x_data.size() == y_data.size() )
    {
        std::string x_data_string = PythonPlot::commandLinePlotData(x_data);
        std::string y_data_string = PythonPlot::commandLinePlotData(y_data);
        PythonPlot plot = PythonPlot(x_data_string,y_data_string,x_label,y_label,legend_data,title_data,save_file_name,x_limits,y_limits);
        plot.plot();
    }
    else
    {
        throw Error::XYDataVectorsNotSameSize;
    }
}

/**
 * Turns a vector<vector<Real>> into a command line set
 * @param data_sets
 * @return 
 */
std::string PythonPlot::commandLinePlotData(const std::vector<std::vector<Real>> &data_sets)
{
    std::string data_string = "";
    size_t number_data_sets = data_sets.size();
    
    //print the data
    for( int index = 0; index < number_data_sets; index++ )
    {
        //grab the data string for the data set
        data_string += PythonPlot::commandLinePlotData(data_sets[index]);
        
        //add a # separator between data sets
        if( index < number_data_sets - 1 )
        {
            data_string += " # ";
        }
        
    }  

    return data_string;
}


/**
 * Turns a vector<Real> into a command line set
 * @param data_set
 * @return 
 */
std::string PythonPlot::commandLinePlotData(const std::vector<Real> &data_set)
{
    std::stringstream data_stream;
    data_stream << std::scientific;
    
    //print the data
    for( std::size_t index = 1; index < data_set.size(); index++ )
    {
        Real data = data_set[index];
        std::string data_str = std::to_string(data) + " ";
        data_stream << data_str;        
    }  

    std::string data_string = data_stream.str();    
    return data_string;
            
}

/**
 * Code to run the command line parameters
 */
void PythonPlot::plot()
{
    
    std::string program = PYTHON_PLOT_SCRIPT;  
    
    
    std::string x_command = " --xdata=\"" + _x_data + "\" ";
    std::string y_command = " --ydata=\"" + _y_data + "\" ";
    std::string axis_labels = " --xlabel='" + _x_label + "' --ylabel='" + _y_label + "' ";
    std::string title = " --title='" + _title + "' ";
    std::string save_file = " --saveplot='" + _save_file + "' ";
    std::string legend = " --legend=\"" + _legend + "\" ";
    
    std::string limits = "";
    
    if(_x_limits.first != _x_limits.second)
    {
        limits += " --xlimits='" + std::to_string(_x_limits.first ) + " " + std::to_string(_x_limits.second ) + "' ";
    }
    
    if(_y_limits.first != _y_limits.second)
    {
        limits += " --ylimits='" + std::to_string(_y_limits.first ) + " " + std::to_string(_y_limits.second ) + "' ";
    }
    
    std::string command_line_plot = program + x_command + y_command  + axis_labels + legend + title + save_file + limits;

    
    this->printToLogFile(command_line_plot);
}

void PythonErrorPlot::plotData(const std::vector<std::tuple<Real,Real,Real> > &data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits )
{
    std::vector<Real> x_data;
    std::vector<Real> y_data;
    std::vector<Real> error_data;
    
    for(int x=0; x < data.size(); x++ )
    {
        x_data.push_back(std::get<0>(data[x]));
        y_data.push_back(std::get<1>(data[x]));
        error_data.push_back(std::get<2>(data[x]));
    }
    PythonErrorPlot::plotData(x_data,y_data,error_data,x_label,y_label,legend_data,title_data,save_file_name, x_limits, y_limits);    
}

void PythonErrorPlot::plotData(const std::vector<Real> &x_data,const std::vector<Real> &y_data,const std::vector<Real> &error_data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits)
{
    
    if( x_data.size() == y_data.size() && x_data.size() == error_data.size() )
    {
        std::string x_data_string = PythonErrorPlot::commandLinePlotData(x_data);
        std::string y_data_string = PythonErrorPlot::commandLinePlotData(y_data);
        std::string error_data_string = PythonErrorPlot::commandLinePlotData(error_data);
        
        PythonErrorPlot plot = PythonErrorPlot(x_data_string,y_data_string,error_data_string, x_label,y_label,legend_data,title_data,save_file_name);
        plot.plot();
    }
    else
    {
        throw Error::XYDataVectorsNotSameSize;
    }  
    
    
}

PythonErrorPlot::PythonErrorPlot(const std::string &x_data, const std::string &y_data,  const std::string &error_data, const std::string &x_label, const std::string &y_label, const std::string &legend_data, const std::string &title_data, const std::string &save_file_name, const std::pair<Real,Real> &x_limits, const std::pair<Real,Real> &y_limits ) :
PythonPlot(x_data, y_data, x_label, y_label, legend_data, title_data, save_file_name, x_limits, y_limits)
{
    _error_data = error_data;
}

void PythonErrorPlot::plot()
{
    std::string program = PYTHON_PLOT_SCRIPT;  
    
    std::string x_command = " --xdata=\"" + _x_data + "\" ";
    std::string y_command = " --ydata=\"" + _y_data + "\" ";
    std::string error_command = " --errordata=\"" + _error_data + "\" ";
    std::string axis_labels = " --xlabel='" + _x_label + "' --ylabel='" + _y_label + "' ";
    std::string title = " --title='" + _title + "' ";
    std::string save_file = " --saveplot='" + _save_file + "' ";
    std::string legend = " --legend=\"" + _legend + "\" ";
    std::string command_line_plot = program + x_command + y_command + error_command  + axis_labels + legend + title + save_file;
    this->printToLogFile(command_line_plot);
}


void PythonPlot::printToLogFile(const std::string &command)
{
    
    
    
    if(PythonPlot::_log_file != "")
    {
        //Save the data
        std::ofstream log_file;
        log_file.open (PythonPlot::_log_file, std::fstream::app);
        log_file << command << std::endl;
        log_file.close(); 
        
        //Create the editable figure caller
        /*std::cout<< this->_save_file + ".sh" << std::endl;
        std::ofstream figure_file;
        figure_file.open ( this->_save_file + ".sh" , std::fstream::out);
        figure_file << "#!/bin/bash" << std::endl;
        figure_file << "$( $( sed  '"<< this->_log_counter << "q;d' " << "graph_log.log )" << " --figure )" << std::endl;
        figure_file.close();
        this->_log_counter++;*/
    }
}

//This function creates all of the plots at one time using the log_file and avoilds the 
//command line long command issues;
void PythonPlot::createPlots()
{
    std::string program = PYTHON_PLOT_SCRIPT;  
    std::string command = program + " inputfile " + PythonPlot::_log_file;
    //exec(command);   
}

std::string PythonPlot::_log_file = "";
int PythonPlot::_log_counter = 1;
