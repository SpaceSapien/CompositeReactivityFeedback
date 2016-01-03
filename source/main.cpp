//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include <math.h>
#include <iomanip>
#include <ctime>
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "ExplicitSolverSettings.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InfiniteCompositeReactor.h"




int main(int argc, char** argv) 
{
   /* int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);
*/
   //   MPI_Finalize();     
//}
    
    InfiniteCompositeReactor reactor = InfiniteCompositeReactor();
    reactor.simulate();

}

std::vector<Real> getPowerDistribution(std::vector<Dimension> radial_points, Real kernel_radius, Real total_power_density)
{
    auto number_points = radial_points.size();
    
    std::vector<Real> power_distribution = std::vector<Real>();
    power_distribution.reserve( number_points);
    Real large_radius = radial_points.back();
    Real kernel_power = total_power_density * pow(large_radius,3)/pow(kernel_radius,3);
    
    for( long index = 0; index < number_points; ++index)
    {
        Dimension radial_point = radial_points[index];
        Real power;
        
        if(radial_point < kernel_radius )
        {
            power = kernel_power;           
        }
        else
        {
            power = 0;           
        }
        
        power_distribution.push_back( power );
    }
    
    return power_distribution;
}



std::string exec(const char* cmd) 
{
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    
    if (!pipe) 
    {
        return "ERROR";
    }
    
    char buffer[128];
    
    std::string result = "";
    
    while (!feof(pipe.get())) 
    {
        if (fgets(buffer, 128, pipe.get()) != NULL)
        {
            result += buffer;
        }
    }
    return result;
}
