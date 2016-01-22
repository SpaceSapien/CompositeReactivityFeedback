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


std::string exec(const std::string command, const bool &print_command,const bool &print_output) 
{
    if(print_command)
    {
        std::cout<<command<<std::endl;
    }
    
    
    const char* cmd = command.c_str();
    
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
            
            if(print_output)
            {
                std::cout<<buffer;
            }
        }
    }
    
    if(print_output)
    {
        std::cout<<std::endl;
    }
    
    
    return result;
}

std::string doubleToScientificString(double value)
{
    std::stringstream output_file_stream;
    output_file_stream << std::scientific << value;
    return output_file_stream.str();
    
}