/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "EnumsAndFunctions.h"
#include <string>
#include <cmath>

std::string getMaterialName(const Materials &material)
{
    switch(material)
    {
        case Materials::U :
        {
            return "U";
        }
        
        case Materials::UO2 :
        {
            return "UO2";
        }
        case Materials::UC :
        {
            return "UC";
        }
        case Materials::UN :
        {
            return "UN";
        }
        case Materials::U3Si :
        {
            return "U3Si";
        }
        case Materials::SiC :
        {
            return "SiC";
        }
        case Materials::C :
        {
            return "C";
        }
        case Materials::Be :
        {
            return "Be";
        }
        case Materials::BeO :
        {
            return "BeO";
        }
        case Materials::ZrB2 :
        {
            return "ZrB2";
        }
        case Materials::W :
        {
            return "W";
        }
        case Materials::B4C :
        {
            return "B4C";
        }
        case Materials::Mo :
        {
            return "Mo";
        }
        case Materials::Nb :
        {
            return "Nb";
        }
        case Materials::Zr :
        {
            return "Zr";
        }
        case Materials::Graphene :
        {
            return "Graphene";
        }
        case Materials::ZrO2 :
        {
            return "ZrO2";
        }
        default :
        {
            return "Unknown Material";
        }
    }
};

Materials getMaterialFromName(const std::string &material_name)
{
    if(material_name == "U")
    {
        return Materials::U;
    }
    if(material_name == "UO2")
    {
        return Materials::UO2;
    }   
    if(material_name == "UN")
    {
        return Materials::UN;
    }
    if(material_name == "UC")
    {
        return Materials::UC;
    }
    if( material_name == "U3Si" )
    {
        return Materials::U3Si;
    }
    if( material_name == "SiC" )
    {
        return Materials::SiC;
    }
    if( material_name == "C" )
    {
        return Materials::C;
    }
    if( material_name == "Be" )
    {
        return Materials::Be;
    }
    if( material_name == "BeO" )
    {
        return Materials::BeO;
    }
    if( material_name == "ZrB2" )
    {
        return Materials::ZrB2;
    }
    if( material_name == "W" )
    {
        return Materials::W;
    }
    if( material_name == "B4C" )
    {
        return Materials::B4C;
    }    
    if( material_name == "Mo" )
    {
        return Materials::Mo;
    }
    if( material_name == "Nb" )
    {
        return Materials::Nb;
    }
    if( material_name == "Zr" )
    {
        return Materials::Zr;
    }
    if( material_name == "ZrO2" )
    {
        return Materials::ZrO2;
    }
    throw -1;
};

