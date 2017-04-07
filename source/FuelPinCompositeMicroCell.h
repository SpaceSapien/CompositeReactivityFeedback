/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FuelPinCompositeMicroCell.h
 * Author: chris
 *
 * Created on April 6, 2017, 6:18 PM
 */

#ifndef FUELPINCOMPOSITEMICROCELL_H
#define FUELPINCOMPOSITEMICROCELL_H
#include "CompositeMicroCell.h"
#include "FuelPinCompositeMicroCell.h"
#include "FuelPinReactor.h"

class FuelPinReactor;

class FuelPinCompositeMicroCell : public CompositeMicroCell
{
public:
    FuelPinCompositeMicroCell(FuelPinReactor* reactor, const Real &temperature);
    virtual ~FuelPinCompositeMicroCell();
private:

};

#endif /* FUELPINCOMPOSITEMICROCELL_H */

