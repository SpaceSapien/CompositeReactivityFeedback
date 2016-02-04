/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputFileParser.cpp
 * Author: chris
 * 
 * Created on January 22, 2016, 6:02 PM
 */
#include "InputDataFunctions.h"
#include "InputFileParser.h"
#include <tuple>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <iostream>

std::vector<std::string> InputFileParser::split(const std::string &str,const std::string &sep)
{
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current=strtok(cstr,sep.c_str());
    
    while(current!=NULL)
    {
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;
}


std::string InputFileParser::getInputFileTextEntry(const std::string &name)
{
    std::string base_command = "cat " + _input_file_name + " | perl -ne ";
    std::string parameter_regex = "'/^" + name + ":[\\s]+?([^#]+)[\\s]+?.*$/ && print $1'";
    std::string command = base_command + parameter_regex;
    std::string value = exec(command);
    
    return value;
}

Real InputFileParser::getInputFileParameter(const std::string &name, const Real &default_value)
{
    std::string input_file_data = this->getInputFileTextEntry(name);
    
    if(input_file_data == "")
    {
        return default_value;
    }
    else
    {
        return std::stod(input_file_data);
    }
}
std::vector<Real> InputFileParser::getInputFileParameter(const std::string &name, const std::vector<Real> &default_value)
{
    std::string input_file_data = this->getInputFileTextEntry(name);
    
    
    if(input_file_data == "")
    {
        return default_value;
    }
    else
    {
        std::vector<std::string> split_string = this->split(input_file_data," ");
        std::vector<Real> return_vector;
        
        for( size_t index = 0; index < split_string.size(); index++)
        {
            Real data_piece = std::stod(split_string[index]);
            return_vector.push_back(data_piece);
            
        }
        
        return return_vector;
    }
}



std::vector<Materials> InputFileParser::getInputFileParameter(const std::string &name, const std::vector<Materials> &default_value)
{
    std::string input_file_data = this->getInputFileTextEntry(name);
    
    
    if(input_file_data == "")
    {
        return default_value;
    }
    else
    {
        std::vector<std::string> split_string = this->split(input_file_data," ");
        std::vector<Materials> return_vector;
        
        for( size_t index = 0; index < split_string.size(); index++)
        {
            Materials material = getMaterialFromName(split_string[index]);
            return_vector.push_back(material);
            
        }
        
        return return_vector;
    }
}
std::string InputFileParser::getInputFileParameter(const std::string &name, const std::string &default_value)
{
    std::string input_file_data = this->getInputFileTextEntry(name);
    
    
    if(input_file_data == "")
    {
        return default_value;
    }
    else
    {
        return input_file_data;
    }
}
bool InputFileParser::getInputFileParameter(const std::string &name, const bool &default_value)
{
    
    std::string input_file_data = this->getInputFileTextEntry(name);
    
    
    if(input_file_data == "")
    {
        return default_value;
    }
    else
    {
        if( input_file_data == "true" )
        {
            return true;
        }
        else if(input_file_data == "false" )
        {
            return false;
        }
        else
        {
            std::cerr<< name << " in the input file " << this->_input_file_name << " has an unknown value assume true.";
            return true;
        }
    }
}

InputFileParser::InputFileParser(std::string input_file_path)
{
    this->_input_file_name = input_file_path;
}

InputFileParser::InputFileParser() {}