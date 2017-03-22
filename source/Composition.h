/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Composition.h
 * Author: chris
 *
 * Created on March 17, 2017, 9:10 PM
 */

#ifndef COMPOSITION_H
#define COMPOSITION_H
#include <vector>
#include <tuple>


template<typename T> 

struct CompositionComponent
{
    CompositionComponent(T object, double amount) : _object(object), _amount(amount) {}
    
    T _object;    
    double _amount;
    double _fraction;
};

template<typename T> 

class Composition 
{
    
public:
    Composition() {}
    virtual ~Composition() {}
    
    void addComponent(const T &object ,const double &amount)
    {
        CompositionComponent<T> new_component(object,amount);        
                
        _data.push_back(new_component);
        
        renormalizeData();
    }
    
    std::size_t size() const
    {
        return _data.size();
    }
    
    CompositionComponent<T> operator[](unsigned int i) const
    { 
        return _data[i];
    }
  
 
    
    
    
private:
    
    void renormalizeData()
    {
        double total;
        
        for( auto iterator = _data.begin(); iterator != _data.end(); ++iterator )
        {
            total += (*iterator)._amount;
        }
        
        for( auto iterator = _data.begin(); iterator != _data.end(); ++iterator )
        {
            iterator->_fraction = iterator->_amount / total;
        }
    }
    
    std::vector<CompositionComponent<T>> _data;

};

#endif /* COMPOSITION_H */

