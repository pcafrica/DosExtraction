/* C++11 */

/**
 * @file   dosModel.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Mathematical model for Density of States extraction.
 *
 */

#ifndef DOSMODEL_H
#define DOSMODEL_H

#include "charge.h"
#include "csvParser.h"
#include "factory.h"
#include "numerics.h"
#include "paramList.h"
#include "quadratureRule.h"
#include "solvers.h"
#include "typedefs.h"

#include "gnuplot-iostream.h"

#include <chrono>    // Timing.
#include <iomanip>    // setf and precision.
#include <limits>    // NaN.

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
        explicit
        DosModel (const ParamList &);
        
        /**
         * @brief Destructor (defaulted).
         */
        virtual
        ~DosModel() = default;
        
        /**
         * @brief Perform the simulation.
         * @param[in] config             : the GetPot configuration object;
         * @param[in] input_experim      : the file containing experimental data;
         * @param[in] output_directory   : directory where to store output files;
         * @param[in] output_plot_subdir : sub-directory where to store @ref Gnuplot files;
         * @param[in] output_filename    : prefix for the output filename.
         */
        void
        simulate (const GetPot &, const std::string &,
                  const std::string &,
                  const std::string &, const std::string &);
                  
        /**
         * @brief Perform post-processing.
         * @param[in]  config           : the GetPot configuration object;
         * @param[in]  output_filename  : prefix for the output filename;
         * @param[in]  input_experim    : the file containing experimental data;
         * @param[out] output_info      : output file containing infos about the simulation;
         * @param[out] output_CV        : output file containing infos about capacitance-voltage data;
         * @param[in]  A_semic          : area of the semiconductor @f$ \left[ m^{-2} \right] @f$;
         * @param[in]  C_sb             : stray capacitance (see @ref ParamList) @f$ \left[ F \right] @f$;
         * @param[in]  x                : the mesh;
         * @param[in]  Dens             : charge-carrier density @f$ \left[ m^{-3} \right] @f$;
         * @param[in]  Phi              : LUMO;
         * @param[in]  semicNodesNo     : number of nodes in the semconductor region;
         * @param[in]  V_simulated      : simulated voltage values @f$ \left[ V \right] @f$;
         * @param[in]  C_simulated      : simulated capacitance values @f$ \left[ F \right] @f$.
         */
        void
        post_process (const GetPot &, const std::string &,
                      const std::string &, std::ostream &, std::ostream &,
                      const Real &, const Real &, const VectorXr &, const MatrixXr &,
                      const MatrixXr &,
                      const Index, const VectorXr &, const VectorXr &);
                      
        /**
         * @brief Save the @ref Gnuplot output files.
         * @param[in] output_directory   : directory where to store output files;
         * @param[in] output_plot_subdir : sub-directory where to store @ref Gnuplot files;
         * @param[in] csv_filename       : .csv file to plot;
         * @param[in] output_filename    : prefix for the output filename.
         */
        void
        save_plot (const std::string &, const std::string &,
                   const std::string &, const std::string &) const;
                   
        /**
         * @brief Defines commands to generate @ref Gnuplot output files.
         * @param[in]  csv_filename : .csv file to plot;
         * @param[out] os           : output stream.
         */
        void
        gnuplot_commands (const std::string &, std::ostream &) const;
        
        /**
         * @name Getter methods
         * @{
         */
        inline const ParamList&
        params() const;
        
        inline const Real&
        error_L2() const;
        
        inline const Real&
        error_H1() const;
        
        inline const Real&
        error_Peak() const;
        
        inline const Real&
        C_acc_experim() const;
        
        inline const Real&
        C_acc_simulated() const;
        
        inline const Real&
        C_dep_experim() const;
        
        /**
         * @}
         */
        
        /**
         * @name Setter methods
         * @{
         */
        inline void
        setSigma (const Real &);
        
        /**
         * @}
         */
        
    private:
        bool initialized_;    /**< @brief bool to determine if @ref DosModel @a param_ has been properly initialized. */
        
        ParamList params_;    /**< @brief The parameter list. */
        
        Real V_shift_;    /**< @brief Peak shift between experimental data and simulated values @f$ [V] @f$. */
        
        Real error_L2_;    /**< @brief @f$ L^2 @f$-distance between experimental and simulated capacitance values. */
        Real error_H1_;    /**< @brief @f$ H^1 @f$-distance between experimental and simulated capacitance values. */
        Real error_Peak_;    /**< @brief Distance between the ordinates in the peak of experimental and simulated derivative of capacitance values with respect to the gate potential. */
        
        Real C_acc_experim_;    /**< @brief Experimental accumulation capacitance, used for automatic fitting @f$ [F] @f$. */
        Real C_acc_simulated_;    /**< @brief Simulated accumulation capacitance, used for automatic fitting @f$ [F] @f$. */
        Real C_dep_experim_;    /**< @brief Experimental depletion capacitance, used for automatic fitting @f$ [F] @f$. */
};

// Implementations.
inline const ParamList&
DosModel::params() const
{
    return params_;
}

inline const Real&
DosModel::error_L2() const
{
    return error_L2_;
}

inline const Real&
DosModel::error_H1() const
{
    return error_H1_;
}

inline const Real&
DosModel::error_Peak() const
{
    return error_Peak_;
}

inline const Real&
DosModel::C_acc_experim() const
{
    return C_acc_experim_;
}

inline const Real&
DosModel::C_acc_simulated() const
{
    return C_acc_simulated_;
}

inline const Real&
DosModel::C_dep_experim() const
{
    return C_dep_experim_;
}

inline void
DosModel::setSigma (const Real & sigma)
{
    assert( sigma >= 0.0 );
    params_.sigma_ = sigma;
}

#endif /* DOSMODEL_H */
