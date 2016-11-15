#include <vector>
#include "RadialMesh.h"
#include <cmath>

RadialMesh::RadialMesh(MicroGeometry* geometry,const int &minimum_nodes, const int &total_nodes, const int &cells_per_zone)
{
    std::size_t number_zones = geometry->_geometry.size();
    
    //We have a minimum number of nodes per zone
    int left_over = total_nodes - number_zones * minimum_nodes;
    //Todo throw error for negative number
    
    //Sphere square adjustment
    Dimension largest_dimension = geometry->_geometry.back().second * std::pow(6.0/M_PI, 1.0/3.0 );;
    int remaining_nodes = left_over;
    
    Dimension zone_start = 0;
    
    for(std::size_t current_zone = 1; current_zone <= number_zones; current_zone++ )
    {
        Dimension zone_end;
        
        //if we're in the last material zone set the zone_end radius to the equivalent sphere volume to the shell
        if(current_zone == number_zones )
        {
            zone_end = largest_dimension;        
        }
        else
        {
            zone_end = geometry->_geometry[current_zone - 1].second;        
        }
        
        Dimension zone_range = zone_end - zone_start;
        Materials current_material = geometry->_geometry[current_zone - 1].first;
        
        int extra_nodes = std::ceil(remaining_nodes * zone_range/ (largest_dimension - zone_start ) );
        remaining_nodes -= extra_nodes;
        int number_nodes_in_zone = extra_nodes + minimum_nodes; 
                
        for( int current_node = 1; current_node <= number_nodes_in_zone; ++current_node)
        {
            Dimension current_node_location = ((static_cast<Dimension>(current_node)-0.5)/static_cast<Dimension>(number_nodes_in_zone))*(zone_range) + zone_start;
            this->pushBack(current_node_location, current_material);
            
            Dimension inner_radius = ((static_cast<Dimension>(current_node)-1)/static_cast<Dimension>(number_nodes_in_zone))*(zone_range) + zone_start;
            Dimension outer_radius = ((static_cast<Dimension>(current_node))/static_cast<Dimension>(number_nodes_in_zone))*(zone_range) + zone_start;
            _inner_radius.push_back(inner_radius);
            _outer_radius.push_back(outer_radius);
            
            Real volume = sphere_volume(outer_radius) - sphere_volume(inner_radius);
            _volume.push_back(volume);
            
            Real inner_surface_area = sphere_surface_area(inner_radius);
            Real outer_surface_area = sphere_surface_area(outer_radius);
            _inner_surface.push_back(inner_surface_area);
            _outer_surface.push_back(outer_surface_area);            
            _zone.push_back(current_zone);
            
            Real fraction_range = ( current_node_location - zone_start) / zone_range ;
            int current_cell = std::floor( cells_per_zone * fraction_range) + 1;
            _cell_in_zone.push_back(current_cell);
        }
        
        zone_start = zone_end;
    }   
   
}


void RadialMesh::pushBack(const Dimension &dimension, const Materials &material)
{
     _materials.push_back( material );
     _position.push_back( dimension );
}

Dimension RadialMesh::getNodeLocation(const int &node) const
{
    return _position[node];
}

Materials RadialMesh::getMaterial(const int &node) const
{
    return _materials[node];
}

std::size_t RadialMesh::numberOfNodes() const
{
    return _position.size();
}

