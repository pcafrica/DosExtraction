/* C++ */

/**
 * @file   quadratureRule.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Quadrature rules.
 *
 */

#ifndef QUADRATURERULE_H
#define QUADRATURERULE_H

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
		friend class GaussianCharge;
		// Now GaussianCharge can access private parameters with no need to copy
		// them through getter methods that slow down the program execution.
		
		/**
		 * @brief Default constructor (deleted since it is required to specify the no. of nodes).
		 */
		QuadratureRule() = delete;
		/**
		 * @brief Constructor.
		 * @param[in] nNodes : the no. of nodes to be used for the quadrature rule.
		 */
		QuadratureRule(const unsigned &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~QuadratureRule() = default;
		
		/**
		 * @brief Apply the quadrature rule in order to compute the nodes and weights.
		 */
		virtual void apply() = 0;
		
		/**
		 * @name Getter methods
		 * @{
		 */
		inline const unsigned & nNodes () const;
		inline const VectorXd & nodes  () const;
		inline const VectorXd & weights() const;
		/**
		 * @}
		 */
		
	protected:
		unsigned nNodes_ ;	/**< @brief The no. of nodes to be used for the quadrature rule. */
		VectorXd nodes_  ;	/**< @brief Vector containing the computed nodes coordinates. */
		VectorXd weights_;	/**< @brief Vector containing the computed weights. */
};

/**
 * @class GaussHermiteRule
 *
 * Compute nodes and weights for the <em>nNodes_</em>-points approximation of
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
		 * @brief Default constructor (deleted since it is required to specify the no. of nodes).
		 */
		GaussHermiteRule() = delete;
		/**
		 * @brief Constructor.
		 * @param[in] nNodes : the no. of nodes to be used for the quadrature rule.
		 */
		GaussHermiteRule(const unsigned &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~GaussHermiteRule() = default;
		
		virtual void apply() override;
		/**
		 * @brief Apply the quadrature rule reading parameters from a configuration file.
		 * @param[in] config : the GetPot configuration object.
		 */
		void apply(const GetPot &);
		
		/**
		 * @brief Compute nodes and weights using an adapted version of the algorithm presented in: @n
		 * William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. @n
		 * Numerical Recipes: The Art of Scientific Computing (3rd edition). @n
		 * Cambridge University Press, New York, NY, USA.
		 */
		void apply_iterative_algorithm(const unsigned & = 1000, const double & = 1.0e-14);
		/**
		 * @brief Compute nodes and weights using an eigendecomposition-based algorithm.
		 */
		void apply_using_eigendecomposition();
};

// Implementations.
inline const unsigned & QuadratureRule::nNodes() const
{
	return nNodes_;
}

inline const VectorXd & QuadratureRule::nodes() const
{
	return nodes_;
}

inline const VectorXd & QuadratureRule::weights() const
{
	return weights_;
}

#endif /* QUADRATURERULE_H */
