/* C++ */

/**
 * @file   numerics.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Generic numeric algorithms.
 *
 */

#ifndef NUMERICS_H
#define NUMERICS_H

#include "typedefs.h"

#include <limits>    // NaN

/**
 * @namespace numerics
 *
 * @brief Namespace for generic numeric algorithms.
 *
 */
namespace numerics
{
    /**
     * @brief Function to sort @ref Eigen vectors.
     * @tparam ScalarType : the scalar type.
     * @param[in] vector  : the vector to be sorted.
     * @returns the sorted vector.
     */
    template<typename ScalarType>
    VectorX<ScalarType> sort(const VectorX<ScalarType> & vector);
    
    /**
     * @brief Function to sort @ref Eigen vectors, keeping track of indexes.
     * @tparam ScalarType : the scalar type.
     * @param[in] vector  : the vector to be sorted.
     * @returns an @ref Eigen vector of pairs: (sorted value, corresponding index in the unsorted vector).
     */
    template<typename ScalarType>
    VectorXpair<ScalarType> sort_pair(const VectorX<ScalarType> & vector);
    
    /**
     * @brief Take the NaNs out from an input vector.
     * @tparam ScalarType : the scalar type.
     * @param[in] vector : the input vector.
     * @returns a vector containing the non-NaN values of @a vector.
     */
    template<typename ScalarType>
    VectorX<ScalarType> nonNaN(const VectorX<ScalarType> & vector);
    
    /**
     * @brief Function to compute approximate integral of @a y with
     * spacing increment specified by @a x, using trapezoidal rule.
     * @param[in] x : the vector of the discrete domain;
     * @param[in] y : the vector of values to integrate.
     * @returns the approximate integral value.
     */
    Real trapz(const VectorXr & x, const VectorXr & y);
    /**
     * @brief Compute the approximate integral of @a y with unit spacing, using trapezoidal rule.
     * @param[in] y : the vector of values to integrate.
     * @returns the approximate integral value.
     */
    Real trapz(const VectorXr & y);
    
    /**
     * @brief Compute the numeric derivative: @f$ \frac{\mathrm{d}y}{\mathrm{d}x} @f$.
     * @param[in] y : the vector of values to differentiate;
     * @param[in] x : the vector of the discrete domain.
     * @returns a vector of the same length as @a y containing the approximate derivative.
     */
    VectorXr deriv(const VectorXr &, const VectorXr &);
    
    /**
     * @brief Linear 1D interpolation. Interpolate @a y, defined at points @a x, at the point @a xNew.
     * @param[in] x    : the vector of the discrete domain;
     * @param[in] y    : the vector of values to interpolate;
     * @param[in] xNew : the point to interpolate at.
     * @returns a scalar containing the interpolated value.
     */
    Real     interp1(const VectorXr &, const VectorXr &, const Real &);
    /**
     * @brief Linear 1D interpolation. Interpolate @a y, defined at points @a x, at the points @a xNew.
     * @param[in] x    : the vector of the discrete domain;
     * @param[in] y    : the vector of values to interpolate;
     * @param[in] xNew : the vector of points to interpolate at.
     * @returns a vector of the same length as @a xNew containing the interpolated values.
     */
    VectorXr interp1(const VectorXr &, const VectorXr &, const VectorXr &);
    
    /**
     * @brief Compute the @f$ L^2 @f$-norm error @b squared between simulated and interpolated experimental values, using @ref trapz.
     * @param[in] interp    : the interpolated values;
     * @param[in] simulated : the simulated values;
     * @param[in] V         : the vector of the electric potential.
     * @returns the value of the @f$ L^2 @f$-norm error.
     */
    Real error_L2(const VectorXr &, const VectorXr &, const VectorXr &);
    
    /**
     * @brief Compute the @f$ L^{\infty} @f$-norm error between simulated and interpolated experimental values.
     * @param[in] interp    : the interpolated values;
     * @param[in] simulated : the simulated values.
     * @returns the value of the @f$ L^{\infty} @f$-norm error.
     */
    Real error_L_inf(const VectorXr &, const VectorXr &);
}

// Implementations.
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
    // Lambda function, specifying a criterion to sort pairs.
    auto sort_pair_lambda = [] (const std::pair<ScalarType, Index> & v1,
                                const std::pair<ScalarType, Index> & v2) -> bool
    {
        return (v1.first < v2.first);
    };
    
    VectorXpair<ScalarType> copy( vector.size() );
    
    for ( Index i = 0; i < copy.size(); ++i )
    {
        copy(i).first = vector(i);
        copy(i).second = i;
    }
    
    std::sort(copy.data(), copy.data() + copy.size(), sort_pair_lambda);
    
    return copy;
}

template<typename ScalarType>
VectorX<ScalarType> numerics::nonNaN(const VectorX<ScalarType> & vector)
{
    // Get number of non-NaN values.
    Index nNonNaN = 0;
    
    for ( Index i = 0; i < vector.size(); ++i )
    {
        if ( !std::isnan(vector(i)) )
        {
            ++nNonNaN;
        }
    }
    
    // Create and fill the new vector.
    VectorX<ScalarType> nonNan = VectorX<ScalarType>::Zero( nNonNaN );
    
    Index k = 0;
    
    for ( Index i = 0; i < vector.size(); ++i )
    {
        if ( !std::isnan(vector(i)) )
        {
            nonNan(k) = vector(i);
            ++k;
        }
    }
    
    return nonNan;
}

#endif /* NUMERICS_H */
