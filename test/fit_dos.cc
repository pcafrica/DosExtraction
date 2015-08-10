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

using namespace constants;

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
        
        // Set number of threads.
        omp_set_num_threads( config("nThreads", (int) nSimulations) );
        
        // Directories names.
        const std::string output_directory   = config("output_directory", "./output" ) + "_fitting/";
        const std::string output_plot_subdir = (std::string) "gnuplot" + "/";
        
        // Fitting parameters.
        Real negative_shift = config("FIT/negative_shift", 1.0) * KB_T;
        Real positive_shift = config("FIT/positive_shift", 1.0) * KB_T;
        
        assert( negative_shift > 0 && positive_shift > 0 );
        
        const Real & nSplits = config("FIT/nSplits", 3);
        
        const unsigned & errorNorm = config("FIT/errorNorm", 2);
        
        // Initialize vectors.
        VectorXr sigma = VectorXr::Zero( 2 * nSplits );
        VectorXr error = VectorXr::Zero( sigma.size() );
        
        VectorXr C_acc_experim   = VectorXr::Zero( sigma.size() );
        VectorXr C_acc_simulated = VectorXr::Zero( sigma.size() );
        VectorXr C_dep_experim   = VectorXr::Zero( sigma.size() );
        
        Real sigmaMin = 0.1 * KB_T;    // Lowest sigma allowed (>=0).
        
        // Create output directories, if they don't exist.
        if ( system( ("exec mkdir " + output_directory + " " + output_directory
                      + output_plot_subdir + " 2> /dev/null").c_str() ) );
                      
        // Create variables to catch error messages inside the parallel region:
        // there are not many ways to throw exceptions outside an OpenMP block.
        std::string ompException;
        bool ompThrewException = false;
        
        // Loop for the simulations.
        
        for ( Index i = 0; i < nSimulations; ++i )
        {
            // Initialize parameter list.
            ParamList params;
            
            if ( nSimulations == parser.nRows() )    // Simulate each row in the input file.
            {
                params = (ParamList) parser.importRow( i + 1 );
            }
            else
            {
                params = (ParamList) parser.importRow( config("indexes", (int) (i + 1), i) );
            }
            
            // Initial guess for sigma (read from the parameter list).
            {
                if ( params.sigma() != sigmaMin )
                {
                    VectorXr temp1 = VectorXr::LinSpaced(nSplits, std::max(params.sigma() - negative_shift, sigmaMin), params.sigma());
                    VectorXr temp2 = VectorXr::LinSpaced(nSplits + 1, params.sigma(), params.sigma() + positive_shift);
                    
                    sigma << temp1, temp2.segment(1, temp2.size() - 1);
                }
                else
                {
                    sigma = VectorXr::LinSpaced(2 * nSplits, params.sigma(), params.sigma() + positive_shift);
                }
            }
            
            // Output filename.
            const std::string output_filename = "output_" + std::to_string( params.simulationNo() );
            
            // Remove possible old files.
            if ( system( ("exec rm -f " + output_directory + output_filename + "* "
                          + output_directory + output_plot_subdir + output_filename + "* 2> /dev/null").c_str() ) );
                          
            // Fitting output file.
            std::ofstream output_fit;
            output_fit.open(output_directory + output_filename + "_fit.txt", std::ios_base::out);
            output_fit.setf(std::ios_base::scientific);
            
            if ( !output_fit.is_open() )
            {
                throw std::ofstream::failure("ERROR: output files cannot be opened or directory does not exist.");
            }
            
            Index iterationsNo = config("FIT/iterationsNo", 3);
            Index minimum = 0;    // The index of the minimum sigma.
            Real sigmaOld = params.sigma();    // Previous minimum value.
            
            // Fitting loop.
            
            for ( Index j = 0; j < iterationsNo; ++j )
            {
                output_fit << "Iteration " << (j + 1) << "/" << iterationsNo << "..." << std::endl;
                
                // Update sigma based on the previous step.
                
                if ( j >= 1 )
                {
                    sigmaOld = sigma(minimum);
                    
                    if ( sigmaOld != sigmaMin )
                    {
                        // Sigma can't be < 0. If so, set it equal to sigmaMin.
                        VectorXr temp1 = VectorXr::LinSpaced(nSplits, std::max(sigmaOld - negative_shift, sigmaMin), sigmaOld);
                        VectorXr temp2 = VectorXr::LinSpaced(nSplits + 1, sigmaOld, sigmaOld + positive_shift);
                        
                        sigma << temp1, temp2.segment(1, temp2.size() - 1);
                    }
                    else
                    {
                        sigma = VectorXr::LinSpaced(2 * nSplits, sigmaOld, sigmaOld + positive_shift);
                    }
                }
                
                // Step 1: find the best sigma.
                #pragma omp parallel for shared(ompException, ompThrewException) private(config) schedule(dynamic, 1)
                
                for ( Index k = 0; k < sigma.size(); ++k )
                {
                    try    // Exception handling inside parallel region.
                    {
                        if ( omp_get_thread_num() == 0 && j == 0 )
                        {
                            if ( i == 0 )
                            {
                                std::cout << std::endl << "Running on " << omp_get_num_threads() << " thread(s)." << std::endl << std::endl;
                            }
                            
                            std::cout << "Performing simulation No. " << params.simulationNo() << " (fitting)..." << std::endl;
                        }
                        
                        // Initialize model.
                        DosModel model;
                        
                        #pragma omp critical
                        {
                            model = (DosModel) params;
                            model.setSigma( sigma(k) );
                        }
                        
                        // Simulate and save output files.
                        model.simulate(config, input_experim, output_directory, output_plot_subdir,
                                       output_filename + "_" + std::to_string(j + 1) + "_" + std::to_string(k + 1));
                                       
                        // Get the desired error.
                        switch ( errorNorm )
                        {
                            case 0:
                                error(k) = model.error_L2();
                                break;
                                
                            case 1:
                                error(k) = model.error_H1();
                                break;
                                
                            case 2:
                                error(k) = model.error_Peak();
                                break;
                        }
                        
                        C_acc_experim(k)   = model.C_acc_experim();
                        C_acc_simulated(k) = model.C_acc_simulated();
                        C_dep_experim(k)   = model.C_dep_experim();
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
                        throw std::runtime_error(ompException);
                    }
                }
                
                // Find the best fitting, i.e. the one with the minimum error.
                error.minCoeff(&minimum);
                
                // Step 2: update C_sb.
                params.setC_sb( params.C_sb() + C_acc_experim(minimum) - C_acc_simulated(minimum) );
                
                // Step 3: update t_semic.
                params.setT_semic( params.eps_semic() * (params.A_semic() / (C_dep_experim(minimum) - params.C_sb())
                                   - params.t_ins() / params.eps_ins()) );
                
                // Print to output.
                output_fit << "\tBest sigma: " << std::setprecision(4) << sigma(minimum) / KB_T;
                output_fit << " (from simulation " << params.simulationNo();
                output_fit << "_" << (j + 1) << "_" << (minimum + 1) << ")" << std::endl;
                
                output_fit.precision(std::numeric_limits<Real>::digits10);
                
                switch ( errorNorm )
                {
                    case 0:
                        output_fit << "\tL2-error: ";
                        break;
                        
                    case 1:
                        output_fit << "\tH1-error: ";
                        break;
                        
                    case 2:
                        output_fit << "\tPeak-error: ";
                        break;
                }
                
                output_fit << error(minimum) << std::endl;
                
                output_fit << "\tC_sb: " << params.C_sb() << std::endl;
                output_fit << "\tt_semic: " << params.t_semic() << std::endl;
                
                if ( j < iterationsNo - 1 )
                {
                    output_fit << std::endl;
                }
                
                // Update fitting parameters.
                if ( sigma(minimum) < sigmaOld )
                {
                    positive_shift = sigmaOld - sigma(minimum);
                }
                else if ( sigma(minimum) > sigmaOld )
                {
                    negative_shift = sigma(minimum) - sigmaOld;
                }
                else if ( sigma(minimum) == sigmaOld )
                {
                    output_fit << "Convergence reached!" << std::endl;
                    
                    break;
                }
            }
            
            output_fit.close();
            
            std::cout << "\t\t\t\tSimulation No. " << params.simulationNo() << " complete!" << std::endl;
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
