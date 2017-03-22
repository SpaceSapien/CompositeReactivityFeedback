//#include <mpi.h>
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
#include "InfiniteHomogenousReactor.h"


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
   
    try
    {
        InputFileParser* input_file_parser;
        
        
        //No Input File
        /*if( argc == 1 )
        {
            //If no input file
            InfiniteCompositeReactor reactor = InfiniteCompositeReactor();
            reactor.simulateTransient();
        }
        //Input File Specified
        else*/ 
        if( argc == 2)
        {
            std::string input_file_name = std::string(argv[1]);

            if( ! file_exists(input_file_name) )
            {
                std::cerr << "input_file_name: " + input_file_name + " Doesn't exist";
                input_file_name = "input/default-input-file.inp";       
                
                if( ! file_exists(input_file_name) )
                {
                    std::cerr << "input_file_name: " + input_file_name + " Doesn't exist";
                    input_file_name = "input/default-input-file.inp";  
                    
                    //Check to make sure the old dire
                    if( ! file_exists(input_file_name) )
                    {
                        std::cerr << "Input File: " + input_file_name + " doesn't exist";
                        throw 1;

                    }

                    //Initialize the input file reader
                    input_file_parser = new InputFileParser( input_file_name );
                    delete input_file_parser;
                }
                
            }

            InfiniteHomogenousReactor reactor(input_file_name);
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
                    
                    //Check to make sure the old dire
                    if( ! file_exists(input_file_name) )
                    {
                        std::cerr << "Input File: " + input_file_name + " doesn't exist";
                        throw 1;

                    }

                    //Initialize the input file reader
                    input_file_parser = new InputFileParser( input_file_name );
                    delete input_file_parser;
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
        /*else if( argc == 4)
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
        }*/
    }
    catch(const std::string &exception)
    {
        std::cerr << "Error: "  << exception << std::endl;
    }
}
