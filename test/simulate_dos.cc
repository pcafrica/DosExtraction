/* C++11 */

/**
 * @file   simulate_dos.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief A test file.
 *
 */

#include "src/dosModel.h"

#include <omp.h>

/**
 *  @brief The @b main function.
 */
int main(const int argc, const char * const * argv, const char * const * envp)
{
    try
    {
        GetPot commandLine(argc, (char **) argv);
        
        const std::string config_directory = commandLine.follow("../config", 2, "-d", "--directory") + "/";
        
        GetPot config = utility::full_path(commandLine.follow("config.pot", 2, "-f", "--file"),
                                           config_directory).c_str();
                                           
        // Input filenames.
        const std::string input_params  = utility::full_path(config("input_params", "input_params.csv" ),
                                          config_directory);
        const std::string input_experim = utility::full_path(config("input_experim", "input_experim.csv" ),
                                          config_directory);
                                          
        CsvParser parser(input_params, config("skipHeaders", true));
        
        // Get number of simulations to be performed.
        Index nSimulations = 0;
        
        {
            bool simulate_all = config("simulate_all", false);
            
            if ( simulate_all == true )
            {
                nSimulations = parser.nRows();
            }
            else
            {
                nSimulations = config.vector_variable_size("indexes");
            }
            
            if ( nSimulations == 0 )
            {
                throw std::ifstream::failure("ERROR: wrong variables \"simulate_all\" and \"indexes\" set in the configuration file.");
            }
        }
        
        Index iterationsNo = config("FIT/iterationsNo", 3);
        
        // Set number of threads.
        omp_set_num_threads( config("nThreads", (int) nSimulations) );
        
        // Directories names.
        const std::string output_directory   = config("output_directory", "./output" ) + "/";
        const std::string output_plot_subdir = (std::string) "gnuplot" + "/";
        
        // Create output directories, if they don't exist.
        if ( system( ("exec mkdir " + output_directory + " " + output_directory
                      + output_plot_subdir + " 2> /dev/null").c_str() ) );
                      
        // Create variables to catch error messages inside the parallel region:
        // there are not many ways to throw exceptions outside an OpenMP block.
        std::string ompException;
        bool ompThrewException = false;
        
        // Loop for the parallel simulations.
        #pragma omp parallel for shared(ompException, ompThrewException) private(config) schedule(dynamic, 1)
        
        for ( Index i = 0; i < nSimulations; ++i )
        {
            try    // Exception handling inside parallel region.
            {
                if ( omp_get_thread_num() == 0 && i == 0 )
                {
                    std::cout << std::endl << "Running on " << omp_get_num_threads() << " thread(s)." << std::endl << std::endl;
                }
                
                // Initialize model.
                ParamList params;
                DosModel model;
                
                #pragma omp critical
                {
                    // Re-initialize configuration file for each thread.
                    config = (GetPot) utility::full_path(commandLine.follow("config.pot", 2, "-f", "--file"),
                                                         config_directory).c_str();
                                                         
                    // Import params.
                    if ( nSimulations == parser.nRows() )    // Simulate each row in the input file.
                    {
                        params = (ParamList) parser.importRow( i + 1 );
                    }
                    else
                    {
                        params = (ParamList) parser.importRow( config("indexes", (int) (i + 1), i) );
                    }
                }
                
                std::stringstream simulationNo;
                simulationNo << std::setw(2) << std::setfill('0') << params.simulationNo();
                
                for ( Index j = 0; j < iterationsNo; ++j )
                {
                        model = (DosModel) params;
                        
                        #pragma omp critical
                        {
                            for ( Index k = 0; k < j; ++k )
                            {
                                std::cout << "  ";
                            }
                            
                            std::cout << "Performing simulation No. " << simulationNo.str() << ", fitting iteration " << (j + 1) << "..." << std::endl;
                        }
                        
                        // Output filename.
                        const std::string output_filename = "output_" + simulationNo.str();
                        
                        // Remove possible old files.
                        if ( system( ("exec rm -f " + output_directory + output_filename + "* "
                                    + output_directory + output_plot_subdir + output_filename + "* 2> /dev/null").c_str() ) );
                                    
                        // Simulate and save output files.
                        model.simulate(config, input_experim, output_directory, output_plot_subdir, output_filename);
                        
                        params.setC_sb( params.C_sb() + model.C_acc_experim() - model.C_acc_simulated() );
                        params.setT_semic( params.eps_semic() * (params.A_semic() / (model.C_dep_experim() - params.C_sb())
                                          - params.t_ins() / params.eps_ins()) );
                }
                
                #pragma omp critical
                std::cout << "\t\t\t\tSimulation No. " << simulationNo.str() << " complete!" << std::endl;
            }
            catch ( const std::exception & genericException )
            {
                #pragma omp critical
                {
                    ompException = genericException.what();
                    ompThrewException = true;
                }
            }
            
            if ( ompThrewException )
            {
                break;
            }
        }
        
        if ( ompThrewException )
        {
            throw std::runtime_error(ompException);
        }
    }
    catch ( const std::exception & genericException )
    {
        std::cerr << genericException.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cout << std::endl;
    utility::print_block("Tasks complete!", std::cout);
    
    return EXIT_SUCCESS;
}
