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
#include "PythonPlot.h"

PythonPlot::PythonPlot()
{
    
}

PythonPlot::PythonPlot(const std::string &x_data, const std::string &y_data, const std::string &x_label, const std::string &y_label,  const std::string &legend_data, const std::string &title_data, const std::string &save_file_name)
{
    _legend = legend_data;
    _save_file = save_file_name;
    _title = title_data;
    _x_data = x_data;
    _x_label = x_label;
    _y_data = y_data;
    _y_label = y_label;
} 

void PythonPlot::plotData(const std::string& x_data, const std::string& y_data, const std::string& x_label, const std::string& y_label, const std::string& legend_data, const std::string& title_data, const std::string& save_file_name)
{
    PythonPlot plot = PythonPlot(x_data,y_data,x_label,y_label,legend_data,title_data,save_file_name);
    plot.plot();
}


void PythonPlot::plotData(const std::vector<std::vector<Real>> &x_data, const std::vector<std::vector<Real>> &y_data, const std::string &x_label, const std::string &y_label,const std::vector<std::string> &legend_data, const std::string &title_data, const std::string &save_file_name)
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
        PythonPlot plot = PythonPlot(x_data_string,y_data_string,x_label,y_label,legend_string,title_data,save_file_name);
        plot.plot();
    }
    else
    {
        throw Error::XYDataVectorsNotSameSize;
    }
}


void PythonPlot::plotData(const std::vector<std::pair<Real,Real> > &data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name)
{
    std::vector<Real> x_data;
    std::vector<Real> y_data;
    
    for(int x=0; x < data.size(); x++ )
    {
        x_data.push_back(data[x].first);
        y_data.push_back(data[x].second);
    }
    PythonPlot::plotData(x_data,y_data,x_label,y_label,legend_data,title_data,save_file_name);
       
}

void PythonPlot::plotData(const std::vector<std::pair<Real, std::vector<Real> > > &data, const std::string &x_label, const std::string &y_label,const std::vector<std::string> &legend_data, const std::string &title_data, const std::string &save_file_name)
{
 
    std::vector< std::vector<Real> > y_data;
    std::vector< std::vector<Real> >x_data;
    std::vector<Real> temp_x_data;
    std::string legend_string = "";
    
    //If there is no legend data passed in
    if(legend_data.size() == 0)
    {
        //create the legend and the vector<vector<Real>> of times
        for( size_t group_index = 0; group_index < data[0].second.size(); group_index++)
        {
            legend_string += "Group-" + std::to_string(group_index + 1) + " ";
        }
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
    
    
    PythonPlot::plotData(x_data, y_data,x_label,y_label,legend_data,title_data,save_file_name);
}

void PythonPlot::plotData(const std::vector<Real> &x_data, const std::vector<Real> &y_data, const std::string &x_label, const std::string &y_label,const std::string &legend_data, const std::string &title_data, const std::string &save_file_name)
{
    if( x_data.size() == y_data.size() )
    {
        std::string x_data_string = PythonPlot::commandLinePlotData(x_data);
        std::string y_data_string = PythonPlot::commandLinePlotData(y_data);
        PythonPlot plot = PythonPlot(x_data_string,y_data_string,x_label,y_label,legend_data,title_data,save_file_name);
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
    std::string data_string;
    
    //print the data
    for( int index = 0; index < data_set.size(); index++ )
    {
        data_string += std::to_string( data_set[index] ) + " ";
    }  

    return data_string;
}

/**
 * Code to run the command line parameters
 */
void PythonPlot::plot()
{
    std::string program = "/home/chris/PycharmProjects/CommandLinePlot/CommandLinePlot.py ";    
    std::string x_command = " --xdata=\"" + _x_data + "\" ";
    std::string y_command = " --ydata=\"" + _y_data + "\" ";
    std::string axis_labels = " --xlabel='" + _x_label + "' --ylabel='" + _y_label + "' ";
    std::string title = " --title='" + _title + "' ";
    std::string save_file = " --saveplot='" + _save_file + "' ";
    std::string legend = " --legend=\"" + _legend + "\" ";
    std::string command_line_plot = program + x_command + y_command  + axis_labels + legend + title + save_file;

    std::cout<<"\n\n"<<command_line_plot<<"\n\n";
    
    system( command_line_plot.c_str() );
}


