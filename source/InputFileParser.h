/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputFileParser.h
 * Author: chris
 *
 * Created on January 22, 2016, 6:02 PM
 */

#ifndef INPUTFILEPARSER_H
#define INPUTFILEPARSER_H
#include <string>
#include <vector>
#include <tuple>
#include "EnumsAndFunctions.h"


class InputFileParser 
{

public:
    
    std::string _input_file_name;    
    Real getInputFileParameter(const std::string &name, const Real &default_value);
    std::vector<Real> getInputFileParameter(const std::string &name, const std::vector<Real> &default_value);
    std::string getInputFileParameter(const std::string &name, const std::string &default_value);
    bool getInputFileParameter(const std::string &name, const bool &default_value);
    std::vector<Materials> getInputFileParameter(const std::string &name, const std::vector<Materials> &default_value);
    int getInputFileParameter(const std::string &name, const int &default_value);
    std::string getInputFileTextEntry(const std::string &name);
    
    InputFileParser(std::string input_file_path);
    InputFileParser();
    
    
private:
    
    std::vector<std::string> split(const std::string &str,const std::string &sep);

};

#endif /* INPUTFILEPARSER_H */

