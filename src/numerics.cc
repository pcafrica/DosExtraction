/* C++11 */

/**
 * @file   numerics.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 */

#include "numerics.h"

Real numerics::trapz(const VectorXr & x, const VectorXr & y)
{
	assert( x.size() == y.size() );
	
	Real integral = 0.0;
	
	#pragma omp parallel for default(shared) reduction(+: integral)
	
	for ( Index i = 0; i < x.size() - 1; ++i ) {
		integral += 0.5 * ( x(i + 1) - x(i) ) * ( y(i) + y(i + 1) );
	}
	
	return integral;
}

Real numerics::trapz(const VectorXr & y)
{
	return trapz(VectorXr::LinSpaced(y.size(), 1, y.size()), y);
}

VectorXr numerics::deriv(const VectorXr & y, const VectorXr & x)
{
	assert( y.size() == x.size() );
	
	Index n = x.size();
	
	VectorXr dy_dx(n);
	
	dy_dx(0) = ( y(1) - y(0) ) / ( x(1) - x(0) );	// Forward difference.
	
	for ( Index i = 1; i < n - 1; ++i ) {
		dy_dx(i) = ( y(i + 1) - y(i - 1) ) / ( x(i + 1) - x(i - 1) );	// Central difference.
	}
	
	dy_dx(n - 1) = ( y(n - 1) - y(n - 2) ) / ( x(n - 1) - x(n - 2) );	// Backward difference.
	
	return dy_dx;
}

Real numerics::interp1(const VectorXr & x, const VectorXr & y, const Real & xNew)
{
	assert( x.size() == y.size() );
	assert( x == sort(x) );	// Initial grid must be monotonic.
	
	if ( xNew < x.minCoeff() || xNew > x.maxCoeff() ) {	// New point cannot be external to the initial grid.
		return std::numeric_limits<Real>::quiet_NaN();
	}
	
	// Index of the last element in "x" lower than xNew.s: "x(index + 1)" is greater or equal than "xNew".
	Index index = std::lower_bound(x.data(), x.data() + x.size(), xNew) - x.data() - 1;
	
	if ( x(index + 1) == xNew ) {	// If "xNew" belongs to the original grid.
		return y(index + 1);
	}
	
	return (Real) ( (xNew - x(index + 1)) * y(index) - (xNew - x(index)) * y(index + 1) ) / (x(index) - x(index + 1));
}

VectorXr numerics::interp1(const VectorXr & x, const VectorXr & y, const VectorXr & xNew)
{
	assert( x.size() == y.size() );
	assert( x == sort(x) );	// Initial grid must be monotonic.
	
	VectorXr yNew = VectorXr::Zero( xNew.size() );
	
	for ( Index i = 0; i < yNew.size(); ++i ) {
		yNew(i) = interp1(x, y, xNew(i));
	}
	
	return yNew;
}

Real numerics::error_L2(const VectorXr & interp, const VectorXr & simulated,
                        const VectorXr & V, const Real & V_shift)
{
	// Number of not-NaN values.
	Index nNotNaN = std::min( (interp   .array() != std::numeric_limits<Real>::quiet_NaN()).count(),
	                          (simulated.array() != std::numeric_limits<Real>::quiet_NaN()).count()
	                        );
	                        
	VectorXr         V_centered = VectorXr::Zero( nNotNaN );
	VectorXr simulated_centered = VectorXr::Zero( nNotNaN );
	VectorXr    interp_centered = VectorXr::Zero( nNotNaN );
	
	{
		Index k = 0;
		
		for ( Index i = 0; i < nNotNaN; ++i ) {
			if ( !std::isnan(interp(i)) && !std::isnan(simulated(i)) ) {
				V_centered(k) = V(i);
				simulated_centered(k) = simulated(i);
				interp_centered(k) = interp(i);
				++k;
			}
		}
	}
	
	return trapz(V_centered.array() - V_shift, (interp_centered - simulated_centered).array().square().matrix());
}
