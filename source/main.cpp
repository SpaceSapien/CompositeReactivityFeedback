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
#include <sstream>
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
    
    //No Input File
    if( argc == 1 )
    {
        //char[] input_file = argv[1]
        std::string input_file_name = "input/default-input-file.inp";
        InfiniteCompositeReactor reactor = InfiniteCompositeReactor(input_file_name);
        reactor.simulate();
    }
    //Input File Specified
    else if( argc == 2)
    {
        std::string input_file_name = std::string(argv[1]);
        
        if( ! file_exists(input_file_name) )
        {
            std::cerr << "input_file_name: " + input_file_name + " Doesn't exist";
            input_file_name = "input/default-input-file.inp";
            
        }
        
        InfiniteCompositeReactor reactor = InfiniteCompositeReactor(input_file_name);
        reactor.simulate();
    }

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

bool file_exists (const std::string &name) 
{
    if (FILE *file = fopen(name.c_str(), "r")) 
    {
        fclose(file);
        return true;
    } 
    else 
    {
        return false;
    }   
}


std::string doubleToScientificString(double value)
{
    std::stringstream output_file_stream;
    output_file_stream << std::scientific << value;
    return output_file_stream.str();
    
}