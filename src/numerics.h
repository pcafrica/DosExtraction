/* C++ */

/**
 * @file   numerics.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Generic numeric algorithms.
 *
 */

#ifndef NUMERICS_H
#define NUMERICS_H

#include "typedefs.h"

#include <limits>	// NaN

/**
 * @namespace numerics
 *
 * @brief Namespace for generic numeric algorithms.
 *
 */
namespace numerics
{
	/**
	 * @brief Function to sort Eigen vectors.
	 * @tparam ScalarType : the scalar type.
	 * @param[in] vector  : the vector to be sorted.
	 * @returns the sorted vector.
	 */
	template<typename ScalarType>
	VectorX<ScalarType> sort(const VectorX<ScalarType> & vector);
	
	/**
	 * @brief Function to sort Eigen vectors, keeping track of indexes.
	 * @tparam ScalarType : the scalar type.
	 * @param[in] vector  : the vector to be sorted.
	 * @returns an Eigen vector of pairs: (sorted value, corresponding index in the unsorted vector).
	 */
	template<typename ScalarType>
	VectorXpair<ScalarType> sort_pair(const VectorX<ScalarType> & vector);
	
	/**
	 * @brief Function to compute approximate integral of @a y with
	 * spacing increment specified by @a x, using trapezoidal rule.
	 * @param[in] x : the vector of the discrete domain;
	 * @param[in] y : the vector of values to integrate.
	 * @returns the approximate integral value.
	 */
	double trapz(const VectorXd & x, const VectorXd & y);
	/**
	 * @brief Compute the approximate integral of @a y with unit spacing, using trapezoidal rule.
	 * @param[in] y : the vector of values to integrate.
	 * @returns the approximate integral value.
	 */
	double trapz(const VectorXd & y);
	
	/**
	 * @brief Compute the numeric derivative: @f$ \frac{\mathrm{d}y}{\mathrm{d}x} @f$.
	 * @param[in] y : the vector of values to differentiate;
	 * @param[in] x : the vector of the discrete domain.
	 * @returns a vector of the same length as @a y containing the approximate derivative.
	 */
	VectorXd deriv(const VectorXd &, const VectorXd &);
	
	/**
	 * @brief Linear 1D interpolation. Interpolate @a y, defined at points @a x, at the point @a xNew.
	 * @param[in] y    : the vector of values to interpolate;
	 * @param[in] x    : the vector of the discrete domain;
	 * @param[in] xNew : the point to interpolate at.
	 * @returns a scalar containing the interpolated value.
	 */
	double   interp1(const VectorXd &, const VectorXd &, const double &);
	/**
	 * @brief Linear 1D interpolation. Interpolate @a y, defined at points @a x, at the points @a xNew.
	 * @param[in] y    : the vector of values to interpolate;
	 * @param[in] x    : the vector of the discrete domain;
	 * @param[in] xNew : the vector of points to interpolate at.
	 * @returns a vector of the same length as @a xNew containing the interpolated values.
	 */
	VectorXd interp1(const VectorXd &, const VectorXd &, const VectorXd &);
	
	/**
	 * @brief Compute the @f$ L^2 @f$-norm error between simulated and interpolated values, using @a trapz.
	 * @param[in] interp    : the interpolated values;
	 * @param[in] simulated : the simulated values;
	 * @param[in] V         : the vector of the electric potential;
	 * @param[in] V_shift   : shift to the electric potential.
	 * @returns the value of the @f$ L^2 @f$-norm error.
	 */
	double error_L2(const VectorXd &, const VectorXd &, const VectorXd &, const double &);
}

template<typename ScalarType>
VectorX<ScalarType> numerics::sort(const VectorX<ScalarType> & vector)
{
	VectorX<ScalarType> copy( vector );
	std::sort(copy.data(), copy.data() + copy.size());
	
	return copy;
}

template<typename ScalarType>
VectorXpair<ScalarType> numerics::sort_pair(const VectorX<ScalarType> & vector)
{
	VectorXpair<ScalarType> copy( vector.size() );
	
	for ( int i = 0; i < copy.size(); ++i ) {
		copy(i).first = vector(i);
		copy(i).second = i;
	}
	
	std::sort(copy.data(), copy.data() + copy.size(),
	[&](std::pair<ScalarType, unsigned> l, std::pair<ScalarType, unsigned> r) {
		return l.first < r.first;
	}
	         );
	         
	return copy;
}

#endif /* NUMERICS_H */
