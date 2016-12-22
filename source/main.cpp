//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <fstream>
#include "MicroSolution.h"
#include "MaterialLibrary.h"
#include "EnumsAndFunctions.h"
#include "MicroGeometry.h"
#include "MicroCell.h"
#include "InputDataFunctions.h"
#include "ReactorKinetics.h"
#include "PythonPlot.h"
#include "InfiniteCompositeReactor.h"
#include "Tally.h"
#include "TallyGroup.h"


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
        //If no input file
        InfiniteCompositeReactor reactor = InfiniteCompositeReactor();
        reactor.simulateTransient();
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
        reactor.simulateTransient();
    }
    // Only do a worth calculation
    // (1) ./executable (2) "worth" (3) input_file
    else if(argc == 3)
    {
        std::string command = std::string(argv[1]);
        
        if(command == "worth")
        {
            std::string input_file_name = std::string(argv[2]);    
            if( ! file_exists(input_file_name) )
            {
                std::cerr << "input_file_name: " + input_file_name + " Doesn't exist";
                input_file_name = "input/default-input-file.inp";            
            }
            
            InfiniteCompositeReactor reactor = InfiniteCompositeReactor(input_file_name);
            reactor.worthStudy();
        }
        else
        {
            std::cerr << "Unknown Command " << command << std::endl;
        }
        
    }
    //We are resuming a run
    // (1) ./executable (2) "resume" (3) float new end time (4) old_results_folder 
    else if( argc == 4)
    {
        std::string command = std::string(argv[1]);
        
        if(command == "resume")
        {
            Real new_end_time = static_cast<Real>(std::stod(std::string(argv[2])));
                        
            std::string old_results_folder = std::string(argv[3]);
            
            if( ! file_exists(old_results_folder) )
            {
                std::cerr << "input_file_name: " + old_results_folder + " Doesn't exist";
            }

            InfiniteCompositeReactor reactor = InfiniteCompositeReactor(old_results_folder, new_end_time);
            reactor.simulateTransient();
        }        
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

std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last-first+1));
}

std::string get_file_text(const std::string& file_path)
{
    std::stringstream file_data;
    std::ifstream infile(file_path);
    std::string line;
    
    while( std::getline(infile,line))
    {
        file_data << line << std::endl;
    }
    
    
    return file_data.str();
}


Real sphere_volume(const Real &radius)
{
    return M_PI * (4.0/3.0) * std::pow(radius,3);
}

Real sphere_surface_area(const Real &radius)
{
    return M_PI * (4.0) * radius*radius;
}

/**
 * 
 * @param vector1  vector of values that should be the same size as the second vector
 * @param vector2
 * @param max_relative_residual the normalized max residual of any element
 * @param average_residual      the avg normalized residual  for all elements
 */
void vector_residuals(const std::vector<Real> &vector1, const std::vector<Real> &vector2, Real &max_relative_residual,Real &average_residual)
{
    max_relative_residual = 0;
    average_residual = 0;
    std::size_t size = vector1.size();
    
    if( size == vector2.size() )
    {
        for(std::size_t index = 0; index < size; ++index)
        {
            Real current_residual;
            
            if(vector1[index] != 0)
            {
                current_residual = std::abs(vector2[index]/vector1[index] - static_cast<Real>(1.0) );                
            }
            else if( vector2[index] == 0 )
            {
                current_residual = 0;
            }
            else
            {
                //might throw an exception here... maybe
                current_residual = 0;
            }
            
            average_residual += current_residual;
            
            if(current_residual > max_relative_residual)
            {
                max_relative_residual = current_residual;
            }
        }
        
        average_residual /= size;
        
    }
    else
    {
        throw std::string("Vectors are different Sizes ") + std::to_string(__LINE__) + " " + __FILE__;
    }
}

Real vector_max(const std::vector<Real> &vector)
{
    Real max = vector[0];

    for(std::size_t index = 1; index < vector.size(); ++index )
    {
        if( vector[index] > max)
        {
            max = vector[index];
        }
    }
    
    return max;    
}