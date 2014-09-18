/* C++ */

/**
 * @file   dosModel.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Mathematical model for Density of States extraction.
 *
 */

#ifndef DOSMODEL_H
#define DOSMODEL_H

#include "charge.h"
#include "csvParser.h"
#include "numerics.h"
#include "paramList.h"
#include "quadratureRule.h"
#include "solvers.h"
#include "typedefs.h"

#include "gnuplot-iostream.h"

#include <chrono>	// Timing
#include <limits>	// NaN

/**
 * @class DosModel
 *
 * @brief Class providing methods to process a simulation to extract the Density of States
 * starting from a parameter list.
 *
 */
class DosModel
{
	public:
		/**
		 * @brief Default constructor.
		 */
		DosModel();
		
		/**
		 * @brief Explicit conversion constructor.
		 * @param[in] params : a parameter list.
		 */
		explicit DosModel(const ParamList &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~DosModel() = default;
		
		/**
		 * @brief Getter method.
		 */
		inline const ParamList & params() const;
		
		/**
		 * @brief Perform the simulation.
		 * @param[in] config             : the GetPot configuration object;
		 * @param[in] input_experim      : the file containing experimental data;
		 * @param[in] output_directory   : directory where to store output files;
		 * @param[in] output_plot_subdir : sub-directory where to store @ref Gnuplot files;
		 * @param[in] output_filename    : prefix for the output filename.
		 */
		void simulate(const GetPot &, const std::string &, const std::string &,
		              const std::string &, const std::string &) const;
		              
		/**
		 * @brief Perform post-processing.
		 * @param[in]  config         : the GetPot configuration object;
		 * @param[in]  input_experim  : the file containing experimental data;
		 * @param[out] output_fitting : output file containing infos about fitting experimental data;
		 * @param[out] output_CV      : output file containing infos about capacitance-voltage data;
		 * @param[in]  x_semic        : the mesh corresponding to the semiconductor domain;
		 * @param[in]  dens           : charge density;
		 * @param[in]  V_simulated    : simulated voltage values;
		 * @param[in]  C_simulated    : simulated capacitance values.
		 */
		void post_process(const GetPot &, const std::string &, std::ostream &, std::ostream &,
		                  const VectorXr &, const VectorXr &, const VectorXr &, const VectorXr &) const;
		                  
		/**
		 * @brief Defines commands to generate @ref Gnuplot output files.
		 * @param[in]  output_CV_filename : output CV filename;
		 * @param[out] os                 : output stream.
		 */
		void gnuplot_commands(const std::string &, std::ostream &) const;
		/**
		 * @brief Save the @ref Gnuplot output files.
		 * @param[in] output_directory   : directory where to store output files;
		 * @param[in] output_plot_subdir : sub-directory where to store @ref Gnuplot files;
		 * @param[in] output_CV_filename : output CV filename;
		 * @param[in] output_filename    : prefix for the output filename.
		 */
		void save_plot(const std::string &, const std::string &, const std::string &, const std::string &) const;
		
	private:
		bool initialized_;	/**< @brief bool to determine if @ref DosModel @a param_ has been properly initialized. */
		
		ParamList params_;	/**< @brief The parameter list. */
};

// Implementations.
inline const ParamList & DosModel::params() const
{
	return params_;
}

#endif /* DOSMODEL_H */
