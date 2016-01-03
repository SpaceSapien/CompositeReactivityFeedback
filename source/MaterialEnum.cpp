/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "EnumsAndFunctions.h"
#include <string>

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
        default :
        {
            return "Unknown Material";
        }
    }
};