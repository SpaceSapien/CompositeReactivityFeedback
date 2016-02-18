/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PythonPlot.h
 * Author: chris
 *
 * Created on December 27, 2015, 7:14 PM
 */

#ifndef PYTHONPLOT_H
#define PYTHONPLOT_H
#include <string>
#include <vector>
#include <tuple>
#include "EnumsAndFunctions.h"

    
#define PYTHON_PLOT_SCRIPT "python python/plot.py " 
    


class PythonPlot 
{
public:
    
    enum Error
    {
        XYDataVectorsNotSameSize
    };
    
    std::string _x_data;
    std::string _y_data;
    
    std::string _x_label;
    std::string _y_label;
    std::string _legend;
    std::string _title;
    std::string _save_file;
    static std::string _log_file; 
    static int _log_counter;
    
    std::pair<Real, Real> _x_limits;
    std::pair<Real, Real> _y_limits;
    
    PythonPlot();
    PythonPlot(const std::string &x_data, const std::string &y_data,                                                 const std::string &x_label, const std::string &y_label,       const std::string &legend_data,                 const std::string &title_data,    const std::string &save_file_name,    const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    void static plotData(const std::string &x_data, const std::string &y_data,                                       const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    void static plotData(const std::vector<std::vector<Real>> &x_data, const std::vector<std::vector<Real>> &y_data, const std::string &x_label="", const std::string &y_label="", const std::vector<std::string> &legend_data={}, const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    void static plotData(const std::vector<Real> &x_data, const std::vector<Real> &y_data,                           const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0}); 
    void static plotData(const std::vector<std::pair<Real,Real> > &data,                                             const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    void static plotData(const std::vector<std::pair<Real, std::vector<Real> > > &data,                              const std::string &x_label="", const std::string &y_label="", const std::vector<std::string> &legend_data={}, const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});

    void static createPlots();
    virtual void plot();
    
protected:
    
    std::string static commandLinePlotData(const std::vector<std::vector<Real>> &data_set);
    std::string static commandLinePlotData(const std::vector<Real> &data_set);
    void printToLogFile(const std::string &command);

};


class PythonErrorPlot : PythonPlot
{
public:
    
    std::string _error_data;
    
    PythonErrorPlot(const std::string &x_data, const std::string &y_data,  const std::string &error_data, const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});   
    
    void static plotData(const std::vector<std::tuple<Real,Real,Real> > &data, const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="", const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    void static plotData(const std::vector<Real> &x_data, const std::vector<Real> &y_data, const std::vector<Real> &error_data, const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="", const std::string &title_data="", const std::string &save_file_name="", const std::pair<Real,Real> &y_dimensions = {0,0}, const std::pair<Real,Real> &x_dimensions = {0,0});
    
    virtual void plot();
};

#endif /* PYTHONPLOT_H */

