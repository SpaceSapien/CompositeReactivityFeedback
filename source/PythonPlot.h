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
#include "EnumsAndFunctions.h"

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
    
    PythonPlot();
    PythonPlot(const std::string &x_data, const std::string &y_data,                                                 const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="");   
    void static plotData(const std::string &x_data, const std::string &y_data,                                       const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="");
    void static plotData(const std::vector<std::vector<Real>> &x_data, const std::vector<std::vector<Real>> &y_data, const std::string &x_label="", const std::string &y_label="", const std::vector<std::string> &legend_data={}, const std::string &title_data="", const std::string &save_file_name="");
    void static plotData(const std::vector<Real> &x_data, const std::vector<Real> &y_data,                           const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name=""); 
    void static plotData(const std::vector<std::pair<Real,Real> > &data,                                             const std::string &x_label="", const std::string &y_label="", const std::string &legend_data="",              const std::string &title_data="", const std::string &save_file_name="");
    void static plotData(const std::vector<std::pair<Real, std::vector<Real> > > &data,                              const std::string &x_label="", const std::string &y_label="", const std::vector<std::string> &legend_data={}, const std::string &title_data="", const std::string &save_file_name="");

    void plot();
    
private:
    
    std::string static commandLinePlotData(const std::vector<std::vector<Real>> &data_set);
    std::string static commandLinePlotData(const std::vector<Real> &data_set);

};

#endif /* PYTHONPLOT_H */

