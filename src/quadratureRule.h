/* C++11 */

/**
 * @file   quadratureRule.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Quadrature rules.
 *
 */

#ifndef QUADRATURERULE_H
#define QUADRATURERULE_H

#include "sandia_rules.hpp"
#include "typedefs.h"

/**
 * @class QuadratureRule
 *
 * Approximate the integral:
 * @f[ \int_{a}^{b} f(x)~\mathrm{d}x @f]
 * with the finite sum:
 * @f[ \sum_{i=1}^{nNodes\_} w_i \cdot f(x_i) @f]
 * where @f$ \{x_i\}_{i=1}^{nNodes\_} @f$ and @f$ \{w_i\}_{i=1}^{nNodes\_} @f$ are called respectively nodes and weights.
 *
 * @brief Abstract class providing a quadrature rule.
 *
 */
class QuadratureRule
{
    public:
        friend class    GaussianCharge;
        friend class ExponentialCharge;
        // Now GaussianCharge and ExponentialCharge can access private parameters with
        // no need to copy them through getter methods that slow down the program execution.
        
        /**
         * @brief Default constructor (deleted since it is required to specify the number of nodes).
         */
        QuadratureRule() = delete;
        /**
         * @brief Constructor.
         * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
         */
        QuadratureRule(const Index &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~QuadratureRule() = default;
        
        /**
         * @brief Apply the quadrature rule in order to compute the nodes and weights.
         */
        virtual void apply() = 0;
        /**
         * @brief Apply the quadrature rule reading parameters from a configuration file.
         * @param[in] config : the GetPot configuration object.
         */
        virtual void apply(const GetPot & config) = 0;
        
        /**
         * @name Getter methods
         * @{
         */
        inline const Index    & nNodes () const;
        inline const VectorXr & nodes  () const;
        inline const VectorXr & weights() const;
        
        /**
         * @}
         */
        
    protected:
        Index    nNodes_ ;    /**< @brief Number of nodes of quadrature. */
        VectorXr nodes_  ;    /**< @brief Vector containing the computed nodes coordinates. */
        VectorXr weights_;    /**< @brief Vector containing the computed weights. */
};

/**
 * @class GaussHermiteRule
 *
 * Compute nodes and weights for the <em>nNodes_</em>-points approximation of:
 * @f[ \int_{-\infty}^{+\infty} w(x)f(x)~\mathrm{d}x @f]
 * where @f$ w(x) = e^{-x^2} @f$.
 *
 * @brief Class derived from @ref QuadratureRule providing the Gauss-Hermite quadrature rule.
 *
 */
class GaussHermiteRule : public QuadratureRule
{
    public:
        /**
         * @brief Default constructor (deleted since it is required to specify the number of nodes).
         */
        GaussHermiteRule() = delete;
        /**
         * @brief Constructor.
         * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
         */
        GaussHermiteRule(const Index &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~GaussHermiteRule() = default;
        
        virtual void apply() override;
        virtual void apply(const GetPot &) override;
        
        /**
         * @brief Compute nodes and weights using an adapted version of the algorithm presented in: @n
         * William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. @n
         * Numerical Recipes: The Art of Scientific Computing (3rd edition). @n
         * Cambridge University Press, New York, NY, USA.
         */
        void apply_iterative_algorithm(const Index & = 1000, const Real & = 1.0e-14);
        /**
         * @brief Compute nodes and weights using an eigendecomposition-based algorithm.
         */
        void apply_using_eigendecomposition();
};

/**
 * @class GaussLaguerreRule
 *
 * Compute nodes and weights for the <em>nNodes_</em>-points approximation of:
 * @f[ \int_{0}^{+\infty} w(x)f(x)~\mathrm{d}x @f]
 * where @f$ w(x) = e^{-x} @f$.
 *
 * @brief Class derived from @ref QuadratureRule providing the Gauss-Laguerre quadrature rule.
 *
 */
class GaussLaguerreRule : public QuadratureRule
{
    public:
        /**
         * @brief Default constructor (deleted since it is required to specify the number of nodes).
         */
        GaussLaguerreRule() = delete;
        /**
         * @brief Constructor.
         * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
         */
        GaussLaguerreRule(const Index &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~GaussLaguerreRule() = default;
        
        /**
         * @brief Auxiliary function to compute @f$ \log\Gamma\left(x\right) @f$
         * @param[in] x : the point to compute the function at.
         * @returns the natural logarithm of the gamma function evaluated at @a x.
         */
        static Real log_gamma(const Real &);
        
        virtual void apply() override;
        virtual void apply(const GetPot &) override;
        
        /**
         * @brief Compute nodes and weights using an adapted version of the algorithm presented in: @n
         * William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. @n
         * Numerical Recipes: The Art of Scientific Computing (3rd edition). @n
         * Cambridge University Press, New York, NY, USA.
         */
        void apply_iterative_algorithm(const Index & = 1000, const Real & = 1.0e-14);
        /**
         * @brief Compute nodes and weights using an eigendecomposition-based algorithm.
         */
        void apply_using_eigendecomposition();
};

// Implementations.
inline const Index & QuadratureRule::nNodes() const
{
    return nNodes_;
}

inline const VectorXr & QuadratureRule::nodes() const
{
    return nodes_;
}

inline const VectorXr & QuadratureRule::weights() const
{
    return weights_;
}

#endif /* QUADRATURERULE_H */
