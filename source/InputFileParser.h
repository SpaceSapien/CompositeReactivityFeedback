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

class InputFileParser 
{

public:
    
    static std::vector<std::string> InputParameters;
    
    
    InputFileParser();
    InputFileParser(const InputFileParser& orig);
    virtual ~InputFileParser();
private:

};

#endif /* INPUTFILEPARSER_H */

