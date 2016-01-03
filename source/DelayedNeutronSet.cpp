/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DelayedNeutronSet.cpp
 * Author: chris
 * 
 * Created on December 8, 2015, 8:04 PM
 */
#include <string>
#include "DelayedNeutronSet.h"
#include <iostream>

DelayedNeutronSet::DelayedNeutronSet() {}



DelayedNeutronSet::DelayedNeutronSet(const std::vector< std::pair<FissionableIsotope,Real>> &fission_talley)
{
    Real fission_total = 0;
    
    //Todo for mixed fission systems
    for(int index = 0; index < fission_talley.size(); index++)
    {
        fission_total += fission_talley[index].second;
    }
    
    this->_delayed_neutrons_per_fission = 0;
    this->_fast_neutrons_per_fission = 0;
    this->_isotope_name = FissionableIsotope::Mixed;
    this->_delayed_neutrons_per_fission_groups.resize(this->numberDelayedNeutronSets(),0);
    this->_neutrons_per_fission = 0;
    
    for(int index = 0; index < fission_talley.size(); index++)
    {
        DelayedNeutronSet set = DelayedNeutronSet(fission_talley[index].first);
        Real set_fraction = fission_talley[index].first/fission_total;
        this->_delayed_neutrons_per_fission += set._delayed_neutrons_per_fission * set_fraction;
        this->_fast_neutrons_per_fission += set._fast_neutrons_per_fission * set_fraction;
        this->_neutrons_per_fission = set._neutrons_per_fission  * set_fraction;
        
        for(int group_index = 0; group_index < numberDelayedNeutronSets(); group_index++ )
        {
            this->_delayed_neutrons_per_fission_groups[group_index] += set._delayed_neutrons_per_fission_groups[group_index] * set_fraction;
        }        
    }
    
    
}

void DelayedNeutronSet::print()
{
    std::cout<<_isotope_name<<"\n"
             << "Delayed Neutrons Per Fission: "<< _delayed_neutrons_per_fission<< "\n"
             << "Fast Neutrons Per Fission: " << _fast_neutrons_per_fission << "\n"
             << "Neutrons Per Fission: " << _neutrons_per_fission << "\n"
             << "Time Constants ";
    
    for( int index = 0; index < numberDelayedNeutronSets(); index++)
    {
        std::cout<< " " << _beta_time_constants[index];
    }
    
    std::cout<< "\n" << "Decay Constants ";
    
    for( int index = 0; index < numberDelayedNeutronSets(); index++)
    {
        std::cout<< " " << _delayed_neutrons_per_fission_groups[index];
    }
    
    std::cout<<"\n";
}

DelayedNeutronSet DelayedNeutronSet::operator=(const DelayedNeutronSet &original) 
{
    _delayed_neutrons_per_fission = original._delayed_neutrons_per_fission;
    _delayed_neutrons_per_fission_groups = original._delayed_neutrons_per_fission_groups;
    _fast_neutrons_per_fission = original._fast_neutrons_per_fission;
    _isotope_name = original._isotope_name;
    _neutrons_per_fission = original._neutrons_per_fission;    
    return *this;
}


DelayedNeutronSet::DelayedNeutronSet(const FissionableIsotope &isotope_name) 
{
    switch(isotope_name)
    {
        
        case FissionableIsotope::U233 :
        {
            _delayed_neutrons_per_fission_groups = {0.00053, 0.00197, 0.00175, 0.00212, 0.00047, 0.00016 };
            _neutrons_per_fission = 2.5;
            break;
        }    
        case FissionableIsotope::U235 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00060, 0.00364, 0.00349, 0.00628, 0.00179, 0.00070 };
            _neutrons_per_fission = 2.43;
            break;
        }
        case U238 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00049, 0.00540, 0.00681, 0.01526, 0.00836, 0.00488 };
            _neutrons_per_fission = 2.492;
            break;
        }
        case Th232 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00143, 0.00776, 0.00843, 0.02156, 0.00838, 0.00204 };
            _neutrons_per_fission = 2.1047;
            break;
        }
        case Pu239 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00024, 0.00176, 0.00126, 0.00207, 0.00065, 0.00022 };
            _neutrons_per_fission = 2.88;
            break;
        }
        case Pu240 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00028, 0.00237, 0.00162, 0.00314, 0.00106, 0.00039 };
            _neutrons_per_fission = 2.897;
            break;
        }
        case Pu241 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00019, 0.00369, 0.00276, 0.00534, 0.00310, 0.00032 };
            _neutrons_per_fission = 2.9453;
            break;
        }
        case Pu242 :
        {
            _delayed_neutrons_per_fission_groups = { 0.00036, 0.00263, 0.00270, 0.00607, 0.00279, 0.00145 };
            _neutrons_per_fission = 2.893;
            break;
        }
        default :
        {
            throw Error::UnknownIsotope;
        }
                
    }
    
    _delayed_neutrons_per_fission = this->sumVector(_delayed_neutrons_per_fission_groups);
    _fast_neutrons_per_fission = _neutrons_per_fission - _delayed_neutrons_per_fission;
}

Real DelayedNeutronSet::sumVector(std::vector<Real> data)
{
    Real sum = 0;
    
    for(std::vector<Real>::iterator it = data.begin(); it != data.end(); ++it)
    {
        sum += *it;
    }
    
    return sum;
}

Real DelayedNeutronSet::getTotalBeta()
{
    return _delayed_neutrons_per_fission/_neutrons_per_fission;
}

int DelayedNeutronSet::numberDelayedNeutronSets()
{
    return _beta_time_constants.size();
}
