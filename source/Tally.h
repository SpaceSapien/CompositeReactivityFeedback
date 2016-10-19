/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Tally.h
 * Author: chris
 *
 * Created on March 5, 2016, 12:51 AM
 */

#ifndef TALLY_H
#define TALLY_H
#include "EnumsAndFunctions.h"

class Tally 
{
public:
    
    std::string _tally_name;
    Real _time;
    std::vector<Real> _energy_bin_tally_counts;
    std::vector<Real> _energies_bins;
    
    Tally();
    
private:

};

#endif /* TALLY_H */

