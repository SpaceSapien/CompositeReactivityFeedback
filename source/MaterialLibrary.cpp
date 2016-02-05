/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 *
 *  Unless stated otherwise all material properties from 
 *  THERMOPHYSICAL PROPERTIES
 *  OF MATERIALS FOR NUCLEAR ENGINEERING
 *  P. Kirillov
 *  State Science Center of the Russian Federation
/* 
 * File:   Materialibrary.cpp
 * Author: chris
 * 
 * Created on November 18, 2015, 2:02 PM
 */
#include <vector>
#include <complex>
#include <math.h>
#include "MaterialLibrary.h"

//Constructor
MaterialLibrary::MaterialLibrary() 
{
    
}



/**
 * Interpolates data from two vectors based on the x value
 * 
 * @param x_array
 * @param y_array
 * @param x
 * @return 
 */
std::pair<Real,Real> MaterialLibrary::interpolateDataAndTemperatureArraysAndDerivative(const std::vector<Real> &x_array,const std::vector<Real> &y_array,const Real &x)
{
    int size = x_array.size();
    
    Real lower = x_array.front();
    
    if( x < lower )
    {
        return std::pair<Real,Real>(y_array[lower],0);
    }
    
    Real upper = x_array.back();
    
    if( x > upper )
    {
        return std::pair<Real,Real>(y_array[upper],0);
    }
    
    int upper_index = size - 1;
    int lower_index = 0;
    
    int index = round(size/2.0);
    bool less_than_far_point;
    bool more_than_close_point;
    bool interpolation_between_points;
    
    do
    {
        
        less_than_far_point = x <= x_array[index + 1];
        more_than_close_point = x >= x_array[index ];
        interpolation_between_points = less_than_far_point && more_than_close_point;
        
        if( ! more_than_close_point )
        {
            upper_index = index;
            index = upper_index - ceil( (upper_index - lower_index)/2.0 );
        }        
        else if( ! less_than_far_point )
        {
            lower_index = index;
            index = lower_index + ceil( (upper_index-lower_index)/2.0 );
        }        
        
    }while( ! interpolation_between_points );
    
    Real x1 = x_array[index];
    Real x2 = x_array[index+1];
    Real y1 = y_array[index];
    Real y2 = y_array[index+1];
    
    //Linear Interpolation
    Real derivative = (y2 - y1)/(x2 - x1);
    Real value = y1 + (x - x1) * derivative;
    
    return std::pair<Real,Real>(value,derivative);
    
}

//Density - kg/m^3
std::pair<Real,Real> MaterialLibrary::getDensityPair(const Materials &material,const Real &T,const Real &dpa)
{
    Real density = 0;
    Real density_derivative = 0;
    
    switch(material)
    {
        case Materials::U :
        {
        
            if(T>237 && T<=942)
            {
                density = 19.36*10e3 - 1.03347*T;
                density_derivative = -1.03347;
            }
            else if(T<1049)
            {
                density = 19.092*10e3 - .9807*T;
                density_derivative = -0.9807;
            }
            else if(T<1405)
            {
                density = 18.447e3 - .5166*T;
                density_derivative = -0.5166;
            }            
            
            return std::pair<Real,Real>(density,density_derivative);
        
        }
   
        case Materials::UO2 :
        {
            density = 10960;
            // TODO
            //rho =  10960/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);
        }        

        case Materials::UN :
        {
            if(T > 298 && T < 2523)
            {
                density = 14420 - 0.2779*T - 4.897*10e-5*T*T;
                density_derivative = -0.2779 -2 * 4.897 * 10e-5 * T;
            }            
            
            return std::pair<Real,Real>(density,density_derivative);
        }

        case Materials::UC :
        {
            density = 13630 * ( 1 - 3.117 * 10e-5 * (T -273.15) - 3.51 * 10e-9 * (T-273.15) * (T-273.15));
            density_derivative = -0.000956826*T - 3.98711;
            return std::pair<Real,Real>(density,density_derivative);
        }
        
        case Materials::C :
        case Materials::Graphene :    
        {
            //valid up to 4000
            density = 1710.0 - 30.0*(T-298.0)/980.0;
            density_derivative = -30.0/980.0;
            return std::pair<Real,Real>(density,density_derivative);
        }
        
    
        case Materials::Be :
        {
            density = 1869.84 - 0.07168*T - 1.6151*10e-5*T*T;
            density_derivative = -0.07168 - 2*1.6151*10e-5*T;
            return std::pair<Real,Real>(density,density_derivative);
        }         
    
        case Materials::BeO :
        {
            //TODO
            density = 2870;
            //2870/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }
    
         
        case Materials::SiC :
        {
            density = 3210;
            //rho = 3210/(( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::W :
        {
            density = 19250;
            //rho = 19250/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::B4C :
        {
            density = 2520;
            //rho = 2520/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::Mo :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm
           //density = 10200/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            density = 10200;
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::Nb :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm
            density = 8550;
            //rho = 8550/ (( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::Zr :
        {
            density = 6499;
            //rho = 6499/(( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }         
    
        case Materials::ZrB2 :
        {
            density = 6090;
            //rho = 6090/(( MaterialProperties.getIntegratedLinearExpansion(material,300,T,radiation) + 1)^3);
            return std::pair<Real,Real>(density,density_derivative);            
        }

        default :
        {
            throw Errors::MaterialNotDefined;
        }
    }
}

//Thermal Conductivity - W/m-K
std::pair<Real,Real> MaterialLibrary::getThermalConductivityPair(const Materials &material,const Real &T,const Real &dpa)
{
    Real thermal_conductivity = 0;
    Real thermal_conductivity_derivative = 0;
    
    switch(material)
    {
        case Materials::U :
        {
            thermal_conductivity = 22.0 + 0.023 * ( T - 273.0);
            thermal_conductivity_derivative = 0.023;
            return std::pair<Real,Real>(thermal_conductivity,thermal_conductivity_derivative);            
        }
        case Materials::UO2 :
        {
            //radiation is measured in MW days/MTU
            Real radiation = 0;
            Real tau = T/1000.0;
            thermal_conductivity = 100.0 / ( 7.5408 + 17.692*tau + 3.6142*tau*tau ) + 6400.0 * exp(-16.35 / ( tau ) ) / pow(tau,2.5);
            Real tau_right = (T+10)/1000;
            Real thermal_conductivity_right = 100.0 / ( 7.5408 + 17.692*tau_right + 3.6142*tau_right*tau_right ) + 6400.0 * exp(-16.35 / ( tau_right ) ) / pow(tau_right,2.5);
            thermal_conductivity_derivative = ( thermal_conductivity_right - thermal_conductivity ) / 10;
            return std::pair<Real,Real>(thermal_conductivity,thermal_conductivity_derivative);                        
        }
        
        case Materials::UN :
        {
            thermal_conductivity = 1.41*pow(T,0.39);
            thermal_conductivity_derivative = 1.41 * 0.39 * pow(T, .39 -1 );
            return std::pair<Real,Real>(thermal_conductivity,thermal_conductivity_derivative);                        
        }
        
        case Materials::UC :
        {
            const std::vector<Real> T_data = {298,  373,  473,  573,  673,  773,  873,  973,  1073, 1173, 1273, 1373, 1473, 1573, 1673, 1773, 1873, 1973, 2073, 2173, 2273, 2373, 2473, 2573, 2673};
            const std::vector<Real> k_data = {25.3, 24.5, 23.6, 23.1, 23.0, 23.1, 23.6, 24.4, 25.6, 27.0, 28.8, 30.9, 33.4, 36.1, 39.2, 42.6, 46.4, 50.4, 54.8, 59.6, 64.6, 70.0, 75.7, 81.7, 81.7};
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }
        case Materials::C :
        {
            const std::vector<Real> T_data = { 298,   400,  500,  530,  615,  733,  863,  989, 1113, 1315, 1423, 1575, 1717, 1890, 1982, 2155, 2318, 2410, 2503, 2615, 2705, 2810, 2933, 3273, 3523, 3773 };
            const std::vector<Real> k_data = { 100.5, 93.5, 86.5, 83.5, 77.2, 71.3, 66.2, 61,  57.4, 51.3, 48.6, 44.9, 43.2, 38.4, 39,   37.6, 34.5, 34.7, 31,   32.9, 29.2, 29.2, 30.4, 26,   18,   6 };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }
        case Materials::Be :
        {
            thermal_conductivity = 202.5 - 0.1723*T + 5.467*10e-5*T*T;            
            thermal_conductivity_derivative = -0.1723 + 2*5.467*10e-5*T;
            return std::pair<Real,Real>(thermal_conductivity,thermal_conductivity_derivative);   
        }
        case Materials::BeO :
        {
            const std::vector<Real> T_data = {303,   353,   373, 423,    473,  573,   673,   773,   873,   973,   1073,  1173,  1273,  1373,  1473,  1573,  1673,  1773,  1873,  1973,  2073,  2173,  2273,  2373,  2473 };
            const std::vector<Real> k_data = {255,   196.7, 177, 122.8,  89.9, 66.9,  54.1,  47.4,  40.8,  35.7,  30.6,  27,    24,    21,    18.8,  16.5,  15.5,  14,    12.5,  12,    11.9,  11.8,  11.7,  11.6,  11.5 };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);            
        }
        
        case Materials::SiC :
        {
            //http://www.nist.gov/data/PDFfiles/jpcrd529.pdf
            thermal_conductivity = 52000*exp(-1.24*10e-5*(T-273.15))/((T-273.15) + 437);
            thermal_conductivity_derivative = exp(-0.000124*T)*(-6.67014*T - 54884.3)/( (T+163.85)*(T+163.85) );
            return std::pair<Real,Real>(thermal_conductivity,thermal_conductivity_derivative); 
        }
        
        case Materials::W :
        {
            //http://www.nist.gov/data/nsrds/NSRDS-NBS-8.pdf

            const std::vector<Real> T_data = {100,  150, 200,  250,  273,  300,  350, 400,  500,  600,  700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,1700,1800,1900, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400,3600 };
            const std::vector<Real> k_data = {235,  210, 197,  186,  182,  178,  170, 162,  149,  139,  133, 128, 124, 121,  118,  115,  113,  111,  109,  107, 105, 103, 102,  100,  98,   96,   94,   92.5, 91.5, 90.5, 90,  89.5 };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
            
        }
        
        case Materials::B4C :
        {
            //http://physics.unm.edu/kenkre/papers/Art88.pdf  taken as a
            //midpoint between the two sets of B4C
            const std::vector<Real> T_data = { 298, 500, 7000, 900, 1100, 1300, 1500, 1700, 2700 };
            const std::vector<Real> k_data = { 36,  15,  13,   11,  10,   9,    8,    7,    5    };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }
        
        case Materials::Mo :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated
            const std::vector<Real> T_data = {298, 1673, 2893};
            const std::vector<Real> k_data = {138, 100,  60  };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }

        case Materials::Nb :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated really need to find trend
            const std::vector<Real> T_data = { 298, 773, 873, 1473, 2741 };
            const std::vector<Real> k_data = { 50,  60,  62,  76,  76    };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }
            
        case Materials::Zr :
        {
            const std::vector<Real> T_data = { 298,  300,  400,  500,  600,  700,  800,  900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2128};
            const std::vector<Real> k_data = { 21.2, 21.2, 19.6, 19.0, 19.0, 19.3, 19.9, 20.6, 21.5, 22.4, 23.5, 24.6, 25.9, 27.2, 28.5, 30.0, 31.5, 33.0, 34.6, 36.3, 36.7}; 
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }

        case Materials::ZrB2 :
        {        
            //Heat conduction mechanisms in hot pressed ZrB2 and
            //ZrB2-SiC Composites
            //Note last point extrapolated
            const std::vector<Real> T_data = { 25, 200,  400, 600, 800, 1000, 1200, 1500, 3516};
            const std::vector<Real> k_data = { 83, 81.5, 78,  77,  74,  72,   72,   71,   68 }; 
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
        }
                
        case Materials::Graphene :
        {        
            const std::vector<Real> T_data = { 100, 200,  400, 800, 4000};
            const std::vector<Real> k_data = { 200, 1000, 2000,1000, 500}; 
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,k_data,T);
            
        }
        
        default :
        {
            throw Errors::MaterialNotDefined;
        }

    }
}

/**
 * 
 * @param material
 * @param material_card_entry  
 * @param doppler_card         The string text entry for the otfdb card entry. Net needed for every material 
 * @param enrichment_fraction  Only needed for materials containing Uranium. Specifies enricment amount
 */
void MaterialLibrary::getMcnpMaterialCard(const Materials &material, const unsigned int &zone, std::string &material_card_entry, std::string &doppler_card,const Real &enrichment_fraction)
{
    std::stringstream material_cards, db_cards;
    
    
    
    
        
    #ifdef LAPTOP

    std::string U238_cs = "92238.66c";
    std::string U235_cs = "92235.66c";

    #elif  PRACTICE_CLUSTER 

    std::string U238_cs = "92238.80c";
    std::string U235_cs = "92235.80c";    

    #endif
        
    Real U235_fraction = enrichment_fraction;
    Real U238_fraction = (1 - enrichment_fraction);
    
    switch(material)
    {
        case Materials::U :
        {
            material_cards << " m" << zone << "  " << U235_cs << "   " << U235_fraction << std::endl;
            material_cards << "     " << U238_cs << "   " << U238_fraction << std::endl;
            break;
           
        }
        case Materials::UO2 :
        {
            material_cards << " m" << zone << "  8016        2          $UO2" << std::endl;
            material_cards << "     " << U235_cs << "   " << U235_fraction << std::endl;
            material_cards << "     " << U238_cs << "   " << U238_fraction << std::endl;
            material_cards << " mt" << zone << " o2-u.27t           $S(a,b) UO2 @ 1200 K" << std::endl;
            material_cards << "     u-o2.27t" << std::endl;
            break;
                            
        }
        
        case Materials::UN :
        {
            material_cards << " m" << zone << "  7014        1          $UN" << std::endl;
            material_cards << "     " << U235_cs << "   " << U235_fraction << std::endl;
            material_cards << "     " << U238_cs << "   " << U238_fraction << std::endl;
            break;
        }
        
        case Materials::UC :
        {
            material_cards << " m" << zone << "  6000        1          $UC" << std::endl;
            material_cards << "     " << U235_cs << "   " << U235_fraction << std::endl;
            material_cards << "     " << U238_cs << "   " << U238_fraction << std::endl;
            break;
        }
        case Materials::C :
        {
            material_cards << " m" << zone << "  6000    1          $Graphite" << std::endl;
            material_cards << " mt" << zone << " grph.22t           $Graphite S(a,b) treatment @ 500 K" << std::endl;
            break;
        }
        case Materials::Be :
        {
            material_cards << " m" << zone << "  4009    1          $Beryllium" << std::endl;
            material_cards << " mt" << zone << " be.25t             $Be S(a,b) treatment @ 800 K" << std::endl;
            break;
        }
        case Materials::BeO :
        {
            material_cards << " m" << zone << "  4009    1          $Beryllium Oxide" << std::endl;
            material_cards << "     8016        1          " << std::endl;
            material_cards << " mt" << zone << " be-o.25t             $Be S(a,b) treatment @ 800 K" << std::endl;
            material_cards << "     o-be.25t" << "             $Oxygen S(a,b) treatment @ 800 K " << std::endl;
            break;        
        }
        
        case Materials::SiC :
        {
            material_cards << " m" << zone << "  14000    1          $SiC" << std::endl;
            material_cards << "     6000        1          " << std::endl;
            break;    
        }
        
        case Materials::W :
        {
            material_cards << " m" << zone << "  74000    1          $Tungsten" << std::endl;
            break;                
        }
        
        case Materials::B4C :
        {
            material_cards << " m" << zone << "  5011    1          $B4C" << std::endl;
            material_cards << "     6000        4          " << std::endl;
            break;
        }
        
        case Materials::Mo :
        {
            material_cards << " m" << zone << "  42000    1          $Mo" << std::endl;
            break;
        }

        case Materials::Nb :
        {
            material_cards << " m" << zone << "  41093    1          $Nb" << std::endl;
            break;
        }
            
        case Materials::Zr :
        {
            material_cards << " m" << zone << "  40000    1          $Zr" << std::endl;
            break;
        }

        case Materials::ZrB2 :
        {        
            material_cards << " m" << zone << "  42000    1          $ZrB2" << std::endl;
            material_cards << "     5011        2          " << std::endl;
            break;
        }
        
        default :
        {
            throw Errors::MaterialNotDefined;
        }

    }
    material_card_entry = material_cards.str();
}



//Specific Heat - J/kg-K
std::pair<Real,Real> MaterialLibrary::getSpecificHeatPair(const Materials &material,const Real &T,const Real &dpa)
{
    Real specific_heat = 0;
    Real specific_heat_derivative = 0;
    
    switch(material)
    {
        case Materials::U :
        {
            if(T < 942)
            {
                specific_heat = 104.82 + 5.3686e-3*T + 10.1823e-5*T*T;
                specific_heat_derivative = 5.3686e-3 + 2*10.1823e-5*T;
            }
            else if(T < 1049)
            {
                specific_heat = 176.4;                
            }
            else
            {
                specific_heat = 156.8;
            }
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
        case Materials::UO2 :
        {
            Real kg_per_mol = (237.0 + 32.0)/1000.0;

            specific_heat = 52.1743 + 87.951*(T/1000) - 84.2411*(T/1000)*(T/1000) + 31.542*pow( (T/1000), 3 ) - 2.6334*pow( (T/1000), 4 ) + 0.71391 * pow( (T/1000), -2 );
            specific_heat_derivative = 87.951/1000.0 - 2*84.2411*(T/1000.0)/1000.0 + 3*31.542*(T/1000.0)*(T/1000.0)/1000.0 - 4*2.6334*pow( (T/1000), 3 )/1000 - 2*0.71391 * pow( (T/1000), -3 );
            specific_heat = (specific_heat/kg_per_mol);
            specific_heat_derivative = (specific_heat_derivative/kg_per_mol);
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
                
        case Materials::UN :
        {
            Real theta = 365.7;
            specific_heat = pow( 0.2029 * (theta/T), 2) * exp(theta/T)/pow( (exp(theta/T) -1 ), 2 ) + 3.766e-5*T + (1.048e9/ (T*T) )*exp(-18081/T);
            
            Real T_right = T + 10.0;
            Real specific_heat_right = pow( 0.2029 * (theta/T_right), 2) * exp(theta/T_right)/pow( (exp(theta/T_right) -1 ), 2 ) + 3.766e-5*T_right + (1.048e9/ (T_right*T_right) )*exp(-18081/T_right);
          
            
            Real kg_per_mol = (237.0 + 14.01)/1000.0;
            specific_heat = specific_heat*1000/kg_per_mol*.19/.76;
            specific_heat_right = specific_heat_right*1000/kg_per_mol*.19/.76;    
            specific_heat_derivative = ( specific_heat_right - specific_heat ) / 10;
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
                   
        case Materials::UC :
        {
            specific_heat = ( 0.2397 - 5.068e-6*T + 1.7604e-8*T*T - 3488.1/(T*T) ) * 1000;
            specific_heat_derivative = ( - 5.068e-6 + 2*1.7604e-8*T + 2*3488.1/(T*T*T) ) * 1000;
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }   

        case Materials::C :
        case Materials::Graphene :
        {   
            const std::vector<Real> T_data =  { 298,  300,  350,  400,  459,   500, 600,   700,  800,   900,  1000,  1100,  1200,  1300,  1400,  1500,  1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300,  2400,  2500,  2600,  2700,  2800,  2900,  3000,  3100,  3200,  3300,  3400,  3500 };
            const std::vector<Real> Cp_data = { 711,  715,  848,  982, 1108,  1220, 1404, 1542, 1648,  1730,  1795,  1848,  1892,  1929,  1961,  1988,  2013, 2035, 2055, 2073, 2090, 2105, 2120, 2134,  2147,  2160, 2.172,  2183,  2195,  2206,  2216,  2227,  2237,  2247,  2257,  2267 };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);              
        }
               
        case Materials::Be :
        {
            specific_heat = (2.1097 + 0.985*10e-3*T - .381*10e5/(T*T))*1000;
            specific_heat_derivative = (0.985*10e-3 + 2*.381*10e5/(T*T*T))*1000;
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
       
            
        case Materials::BeO :
        {
            const std::vector<Real> T_data =  {373,  473,  573,  673,  773,  873,  973,  1073, 1173, 1273, 1373, 1473, 1573, 1673, 1773, 1873, 1973, 2073, 2173, 2273, 2373, 2473};
            const std::vector<Real> Cp_data = {1229, 1464, 1617, 1731, 1826, 1910, 1988, 2014, 2036, 2057, 2079, 2101, 2122, 2144, 2166, 2187, 2209, 2231, 2253, 2274, 2296, 2318};
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);  
        }       
            
        case Materials::SiC :
        {
            //http://www.nist.gov/data/PDFfiles/jpcrd529.pdf
            specific_heat = 1110 + 0.15*T - 425*exp(-.003*(T-273.15));
            specific_heat_derivative = 0.15 + 2.8933 * exp(-0.003*T);
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }       
            
        case Materials::W :
        {
        
            //http://www.nist.gov/data/PDFfiles/jpcrd263.pdf
            if(T>=298 && T<=3500)
            {
                    Real coeifficients[] = { -0.20869, 23.70345, 5.132062, -1.99922, 0.734168 };
                    Real specific_heat_right = 0;
                    
                    for(int index = -1; index<= 3; ++index )
                    {
                        specific_heat = specific_heat + coeifficients[index+1] * pow((T/1000),index);
                        specific_heat_right = specific_heat + coeifficients[index+1] * pow(((T+10)/1000),index);
                    }
                    
                    Real kg_per_mol = 183.85/1000.0;
                    specific_heat = specific_heat / kg_per_mol;
                    specific_heat_right = specific_heat_right / kg_per_mol;
                    
                    specific_heat_derivative = ( specific_heat_right - specific_heat )/10;
            }   
            
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
       
            
        case Materials::B4C :
        {
            specific_heat = 945;
            specific_heat_derivative = 0;
            return std::pair<Real,Real>(specific_heat,specific_heat_derivative);
        }
       
        case Materials::Mo :
        {

            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated
            const std::vector<Real> T_data = { 298, 1673, 2893 };
            const std::vector<Real> Cp_data ={ 240, 330,  420 };
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);            
        }
       
        case Materials::Nb :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated
            const std::vector<Real> T_data =  {298, 1673, 2893};
            const std::vector<Real> Cp_data = {275, 350, 425};
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);                            
        }
                   
        case Materials::Zr :
        {    
            const std::vector<Real> T_data = { 298, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2128};
            const std::vector<Real> Cp_data ={ 285, 285, 298, 210, 320, 331, 340, 350, 360,  370,  307,  313,  320,  327,  335,  344,  354,  366,  378,  392,  396 }; 
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);            
        }
                   
        case Materials::ZrB2 :
        {
        
            //Heat conduction mechanisms in hot pressed ZrB2 and
            //ZrB2-SiC Composites
            //Note last point extrapolated
            const std::vector<Real> T_data = { 100, 200, 300, 400, 500, 600, 1500, 3516};
            const std::vector<Real> Cp_data ={ 425, 510, 540, 565, 585, 600, 675,  825 }; 
            return this->interpolateDataAndTemperatureArraysAndDerivative(T_data,Cp_data,T);            
        }
        default :
        {
            throw Errors::MaterialNotDefined;
        }
       
    }
    
   
}
        

Real MaterialLibrary::getMeltingPoint(const Materials &material)
{
    
    Real melting_temperature = 0;
        
    switch(material)
    {
        case Materials::U :
        {
            melting_temperature = 1405.3;
            break;
        }
        case Materials::UO2 :
        {
            melting_temperature = 3140;
            break;
        }
        case Materials::UN :
        {
           melting_temperature = 2300;
           break;
        }
        case Materials::UC :
        {
            melting_temperature = 2638;
            break;
        }
        case Materials::C :
        case Materials::Graphene :
        {
            melting_temperature = 4000;
           break;
        }
        case Materials::Be :
        {
            melting_temperature =1560;        
            break;
        }
        case Materials::BeO :
        {
        
            melting_temperature = 2780;
            break;
        }
        case Materials::SiC :
        {
            melting_temperature = 3000;

           break;
        }
        case Materials::W :
        {
            melting_temperature = 3695;

           break;
        }
        case Materials::B4C :
        {
            melting_temperature = 2718;

           break;
        }
        case Materials::Mo :
        {
            melting_temperature = 2893;

           break;
        }
        case Materials::Nb :
        {
            melting_temperature = 2741;
           break;
        }
        case Materials::Zr :
        {
            melting_temperature = 2128;
            break;
        }
        case Materials::ZrB2 :
        {
            //NETZCH Poswer
            melting_temperature = 3516;
            break;
        }
        default :
        {
            throw Errors::MaterialNotDefined;
        }

    }
    
    return melting_temperature;
 }

Real MaterialLibrary::getLinearExpansionCoeifficient(const Materials& material, const Real& T, const Real& dpa)
{
    
    Real alpha = 0;
    
    switch(material)
    {
        case Materials::U :
        {
            const std::vector<Real> T_data = {293,     373,     473,     573,     673,     942,     943,     973,     1000,    1049,    1050,    1073,    1173,    1273,    1383,    1406};
            const std::vector<Real> a_data = {10.2e-6, 10.6e-6, 10.8e-6, 10.9e-6, 10.8e-6, 10.2e-6, 11.6e-6, 11.8e-6, 12.0e-6, 12.4e-6, 13.7e-6, 13.9e-6, 14.5e-6, 15.6e-6, 16.5e-6, 16.5e-6 };
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::UO2 :
        {
            const std::vector<Real> T_data = {298,     373,     473,     573,     673,     773,     873,     973,     1073,    1173,    1273,    1373,    1473,    1573,     1673,     1773,     1873,     1973,     2073,     2173,     2273,    2373,    2473,    2573,    2673,    2773,    2873,    2973,    3073,    3100,    3140};
            const std::vector<Real> a_data = {9.76e-6, 9.76e-6, 9.82e-6, 9.90e-6, 10e-6,   10.1e-6, 10.4e-6, 10.5e-6, 10.7e-6, 11e-6,   11.4e-6, 11.9e-6, 12.4e-6, 13e-6,    13.7e-6,  14.4e-6,  15.2e-6,  16.1e-6,  17.0e-6,  18.1e-6,  19.1e-6, 20.3e-6, 21.5e-6, 22.8e-6, 24.1e-6, 25.5e-6, 27e-6,   28.5e-6, 30.1e-6, 30.9e-6, 30.9e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::UN :
        {
            const std::vector<Real> T_data = {298,    373,     473,     573,     673,     773,     873,     973,     1073,    1173,    1273,    1373,    1473,    1573,    1673,    1773,    1873,    1973,    2073,     2173,     2273,     2300   };
            const std::vector<Real> a_data = {7.5e-6, 7.62e-6, 7.76e-6, 7.90e-6, 8.04e-6, 8.19e-6, 8.29e-6, 8.33e-6, 8.47e-6, 8.75e-6, 8.89e-6, 9.03e-6, 9.17e-6, 9.31e-6, 9.45e-6, 9.59e-6, 9.74e-6, 9.86e-6, 10.02e-6, 10.16e-6, 10.30e-6, 10.34e-6 };
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first; 

            break;
        }
        case Materials::UC :
        {
            const std::vector<Real> T_data = {298,     373,     473,     573,     673,     773,     873,     973,     1073,    1173,    1273,    1373,    1473,    1573,    1673,    1773,    1873,    1973,    2073,    2173,    2273,    2373,    2473,    2573,    2638};
            const std::vector<Real> a_data = {10.1e-6, 10.2e-6, 10.3e-6, 10.4e-6, 10.5e-6, 10.7e-6, 10.8e-6, 10.9e-6, 11.0e-6, 11.1e-6, 11.2e-6, 11.4e-6, 11.5e-6, 11.6e-6, 11.7e-6, 11.8e-6, 11.9e-6, 12.1e-6, 12.2e-6, 12.3e-6, 12.4e-6, 12.5e-6, 12.6e-6, 12.8e-6, 12.9e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;
            
            break;
        }
        case Materials::C :
        case Materials::Graphene :
        {
            //note lots of extrapolation in the last number
            const std::vector<Real> T_data = {293,    373,    473,    573,    673,    773,    873,    973,    1073,    1173,    1273,    4000};
            const std::vector<Real> a_data = {4e-6,   4.2e-6, 5e-6,   5e-6,   5.1e-6, 5.2e-6, 5.3e-6, 5.6e-6, 5.6e-6, 5.7e-6,  5.8e-6,   6e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;
            
            break;
        }
        case Materials::Be :
        {
            //last data point is extrapolation
            const std::vector<Real> T_data = {298,    400,    500,    600,    800,    1000,    1200,    1560};
            const std::vector<Real> a_data = {10e-6,  12e-6,  15e-6,  16e-6,  18e-6,  20e-6,   22e-6,   25e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;
       
            break;
        }
        case Materials::BeO :
        {
            //last data point is extrapolation
            const std::vector<Real> T_data = {298,    373,    473,     573,    673,     773,     873,     973,     1073,    1173,    1273,    1373,    1473,    1573,     1673,     1773,     1873,     1973,     2073,     2173,     2273,    2893,      4393};
            const std::vector<Real> a_data = {5.5e-6, 5.6e-6, 6.05e-6, 6.5e-6, 6.95e-6, 7.37e-6, 7.79e-6, 8.19e-6, 8.57e-6, 8.93e-6, 9.27e-6, 9.58e-6, 9.87e-6, 10.12e-6, 10.35e-6, 10.54e-6, 10.70e-6, 10.81e-6, 10.89e-6, 10.93e-6, 10.92e-6, 10.92e-6, 11e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::SiC :
        {
           //last data point is extrapolation
            const std::vector<Real> T_data = {150,     250,     350,     450,     550,     650,     750,     850,     950,     1050,    1150,    1250,    1350,    1450,    1550,    3000 };
            const std::vector<Real> a_data = {0.50e-6, 1.72e-6, 2.69e-6, 3.45e-6, 4.01e-6, 4.42e-6, 4.69e-6, 4.85e-6, 4.93e-6, 4.95e-6, 4.97e-6, 4.97e-6, 5.01e-6, 5.09e-6, 5.27e-6, 5.27e-6 };
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::W :
        {
                //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
                //point extrapolated
                const std::vector<Real> T_data = {298,    1673,    3695};
                const std::vector<Real> a_data = {4.2e-6, 4.65e-6, 5.4e-6};

                alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

           break;
        }
        case Materials::B4C :
        {
            //http://shodhganga.inflibnet.ac.in/bitstream/10603/11623/10/10_chapter%205.pdf? Similar
            //valid for range 298 K to 1773 K
            const std::vector<Real> T_data = {280,     380,     480,     580,     680,     780,     880,     980,     1080,    1180,    1280,    1380,    1480,    1580,    1680,    1780,    1880,    1980,    2080,    2180,    2280,    2380,    2480,     2580,     2680,     2780,     2880,     2980,     3080,     3180,     3280,     3380,     3480 };
            const std::vector<Real> a_data = {2.31e-6, 2.66e-6, 3.04e-6, 3.40e-6, 3.77e-6, 4.13e-6, 4.50e-6, 4.86e-6, 5.23e-6, 5.59e-6, 5.96e-6, 6.32e-6, 6.69e-6, 7.05e-6, 7.42e-6, 7.78e-6, 8.15e-6, 8.52e-6, 8.88e-6, 9.25e-6, 9.61e-6, 9.98e-6, 10.35e-6, 10.71e-6, 10.11e-6, 11.44e-6, 11.81e-6, 12.18e-6, 12.55e-6, 12.91e-6, 13.28e-6, 13.65e-6, 14.01e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::Mo :
        {
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated
            const std::vector<Real> T_data = {298,    1673,    2893};
            const std::vector<Real> a_data = {5.2e-6, 6.3e-6,  7.4e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::Nb :
        {
         
            //http://www.plansee.com/en/Materials-Molybdenum-402.htm note last
            //point extrapolated
            const std::vector<Real> T_data = {300,    1673,    2893};
            const std::vector<Real> a_data = {7e-6,   9.2e-6,  11.4e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

           break;
        }
        case Materials::Zr :
        {
            const std::vector<Real> T_data ={ 298,      300,      400,      500,     600,     700,     800,     900,     1000,    1100,    1200,    1300,    1400,    1500,    1600,    1700,    1800,    1900,    2000,    2100,    2128 };
            const std::vector<Real> a_data ={ 11.46e-6, 11.42e-6, 10.13e-6, 9.50e-6, 9.19e-6, 9.09e-6, 9.10e-6, 9.20e-6, 9.35e-6, 9.53e-6, 8.30e-6, 8.11e-6, 7.90e-6, 7.66e-6, 7.40e-6, 7.11e-6, 6.80e-6, 6.47e-6, 6.12e-6, 5.76e-6, 5.66e-6}; 
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;

            break;
        }
        case Materials::ZrB2 :
        {
            //NETZCH Poswer
            //Data from NETZCH Ceramics Poster, in reality there is a
            //range here...
            const std::vector<Real> T_data = {298,    3516};
            const std::vector<Real> a_data = {5.2e-6, 6.7e-6};
            alpha = this->interpolateDataAndTemperatureArraysAndDerivative(T_data,a_data,T).first;
            
            break;
        }
        default :
        {
            throw Errors::MaterialNotDefined;
        }
    }
    
    
}


/*
 %%
        function alpha_integral = getIntegratedLinearExpansion(material,T_cold,T_hot,radiation)

            integrating_function = @(T)  MaterialProperties.getLinearExpansion(material,T,radiation);

            alpha_integral = integral(integrating_function,T_cold,T_hot);

        end

      
           
    end
 
 
 */

