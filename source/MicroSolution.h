/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Solution.h
 * Author: chris
 *
 * Created on November 26, 2015, 1:15 PM
 */



#ifndef MICROSOLUTION_H
#define MICROSOLUTION_H
#include <string>
#include <vector>
#include "EnumsAndFunctions.h"

class MicroSolution 
{
public:
    
    const static std::string _x_label;
    const static std::string _y_label;
    
    MicroSolution(const std::vector<Dimension> &dimension, const std::vector<Real> &solution, const Real &time);
    int size();
    void plot(const std::string &file_name = "", const Real &min_temperature=0, const Real &max_temperature=0 );
    void static plotSolutions(const std::vector<MicroSolution> &plot_data,const int &number_plots = 0,const std::string &file_name = "");
    void static saveSolutions(const std::vector<MicroSolution> &plot_data_vector, const std::string &save_folder, const std::string &save_file_name="temperature-data.csv");
    
    enum Errors
    {
        SolutionAndGridSizeDifferent
    };
    
    std::vector<Real> _solution;
    std::vector<Dimension> _grid;
    
    Real _time;
    
    MicroSolution static temperatureDifference(MicroSolution* compare1, MicroSolution* compare2);
    
private:
    
    
};

#endif /* SOLUTION_H */

