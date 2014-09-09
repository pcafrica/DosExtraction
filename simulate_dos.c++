/* C++ */

/**
 * @file   simulate_dos.c++
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Main file.
 *
 */

#include "src/dosModel.h"

/**
 *  @brief main() function.
 */

int main(const int argc, const char * const * argv, const char * const * envp)
{
	try {
		GetPot commandLine(argc, (char **) argv);
		GetPot config     ( commandLine.follow("config.pot", 2, "-f", "--file").c_str() );
		
		// Input filenames.
		const std::string input_params  = config("input_params" , "input_params.csv" );
		const std::string input_experim = config("input_experim", "input_experim.csv");
		
		CsvParser parser(input_params, config("hasHeaders", true));
		
		// Get no. of simulations to be performed.
		unsigned nSimulations = 0;
		
		if ( config("simulate_all", false) ) {
			nSimulations = parser.nRows();
		} else {
			nSimulations = config.vector_variable_size("indexes");
		}
		
		if ( nSimulations < 1 ) {
			throw std::ifstream::failure("ERROR: wrong configuration file...");
		}
		
		// Set no. of threads.
		omp_set_num_threads( config("nThreads", (int) nSimulations) );
		
		// Directories names.
		const std::string output_directory   = config("output_directory", "./output" ) + "/";
		const std::string output_plot_subdir = "gnuplot/";
		
		// Clean up and re-create output directories.
		system( ("exec rm -rf " + output_directory + " 2> /dev/null").c_str() );
		system( ("exec mkdir " + output_directory + " " + output_directory + output_plot_subdir + " 2> /dev/null").c_str() );
		
		// Create variables to catch error messages inside the parallel region:
		// there are not many ways to throw exceptions outside an OpenMP block.
		std::string omp_loopException;
		bool thrownException = false;
		
		// Loop for the parallel simulation.
		#pragma omp parallel for shared(omp_loopException, thrownException) private(config)
		
		for ( unsigned i = 0; i < nSimulations; ++i ) {
			try {	// Exception handling inside parallel region.
				// Initialize model.
				DosModel model;
				
				#pragma omp critical
				{
					// Re-initialize configuration file for each thread.
					config = (GetPot) commandLine.follow("config.pot", 2, "-f", "--file").c_str();
					
					// Import params.
					if ( nSimulations == parser.nRows() ) {	// Simulate each row in the input file.
						model = (DosModel) (ParamList) parser.importRow( i + 1 );
					} else {
						model = (DosModel) (ParamList) parser.importRow( config("indexes", (int) (i + 1), i) );
					}
				}
				
				const std::string output_filename = "output_" + std::to_string(model.params().simulationNo());
				
				// Simulate and save output files.
				model.simulate(config, input_experim, output_directory, output_plot_subdir, output_filename);
			} catch ( const std::exception & genericException ) {
				#pragma omp critical
				{
					omp_loopException = genericException.what();
					thrownException = true;
				}
			}
		}
		
		if ( thrownException ) {
			utility::print_block(omp_loopException.c_str(), std::cerr);
			return 1;
		}
	} catch ( const std::exception & genericException ) {
		utility::print_block(genericException.what(), std::cerr);
		return 1;
	}
	
	utility::print_block("Simulation complete!", std::cout);
	
	return 0;
}
