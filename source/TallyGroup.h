/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TallyGroup.h
 * Author: chris
 *
 * Created on November 10, 2016, 9:30 PM
 */

#ifndef TALLYGROUP_H
#define TALLYGROUP_H
#include <vector>
#include <fstream>
#include "Tally.h"
#include "EnumsAndFunctions.h"

class TallyGroup {
public:
    
    Real _time;
    int _zones;
    int _cells_per_zone;
    std::vector<Tally*> _flux_tallies;
    std::vector<Tally*> _fission_tallies;
    std::vector<Tally*> _absorption_tallies;
    
    TallyGroup(const Real &time,const int &zones,const int &cells_per_zone,const std::vector<Tally*> &flux_tallies,const std::vector<Tally*> &fission_tallies,const std::vector<Tally*> &absorption_tallies);
    
    static TallyGroup* MCNPTallyGroupFactory(const std::string &MCTAL_file,const int &zones,const int &cells_per_zone,const Real &time);
    void print(const bool &print_energies = false);
    int size();
    void outputToFile(const std::string &output_file);
    void outputTallyToFile(const std::string &output_file, Tally* tally, const std::string &location);
    
    virtual ~TallyGroup();
    
private:

};

#endif /* TALLYGROUP_H */

