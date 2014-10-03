/* C++11 */

/**
 * @file   simulate_dos.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief A test file.
 *
 */

#include "src/dosModel.h"

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
                                        
      CsvParser parser(input_params, config("hasHeaders", true));
      
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
      
      // Loop for the parallel simulation.
      #pragma omp parallel for shared(ompException, ompThrewException) private(config) schedule(dynamic, 1)
      
      for ( Index i = 0; i < nSimulations; ++i )
        {
          try      // Exception handling inside parallel region.
            {
              // Initialize model.
              DosModel model;
              
              #pragma omp critical
              {
                // Re-initialize configuration file for each thread.
                config = (GetPot) utility::full_path(commandLine.follow("config.pot", 2, "-f", "--file"),
                                                     config_directory).c_str();
                                                     
                // Import params.
                if ( nSimulations == parser.nRows() )      // Simulate each row in the input file.
                  {
                    model = (DosModel) (ParamList) parser.importRow( i + 1 );
                  }
                else
                  {
                    model = (DosModel) (ParamList) parser.importRow( config("indexes", (int) (i + 1), i) );
                  }
              }
              
              std::cout << "Performing simulation No. " << model.params().simulationNo() << "..." << std::endl;
              
              const std::string output_filename = "output_" + std::to_string(model.params().simulationNo());
              
              // Remove possible old files.
              system( ("exec rm -f " + output_directory + output_filename + "* "
                       + output_directory + output_plot_subdir + output_filename + "* 2> /dev/null").c_str() );
                       
              // Simulate and save output files.
              model.simulate(config, input_experim, output_directory, output_plot_subdir, output_filename);
              
              std::cout << "\tSimulation No. " << model.params().simulationNo() << " complete!" << std::endl;
            }
          catch ( const std::exception & genericException )
            {
              #pragma omp critical
              {
                ompException = genericException.what();
                ompThrewException = true;
              }
            }
        }
        
      if ( ompThrewException )
        {
          throw std::runtime_error(ompException);
        }
    }
  catch ( const std::exception & genericException )
    {
      utility::print_block(genericException.what(), std::cerr);
      return EXIT_FAILURE;
    }
    
  std::cout << std::endl;
  utility::print_block("Tasks complete!", std::cout);
  
  return EXIT_SUCCESS;
}
