#include "numerics.h"

namespace numerics
{
	double trapz(const VectorXd & x, const VectorXd & y)
	{
		assert( x.size() == y.size() );
		
		double integral = 0.0;
		
		#pragma omp parallel for default(shared) reduction(+: integral)
		
		for ( int i = 0; i < x.size() - 1; ++i ) {
			integral += 0.5 * ( x(i + 1) - x(i) ) * ( y(i) + y(i + 1) );
		}
		
		return integral;
	}
	
	double trapz(const VectorXd & y)
	{
		return trapz(VectorXd::LinSpaced(y.size(), 1, y.size()), y);
	}
	
	VectorXd deriv(const VectorXd & y, const VectorXd & x)
	{
		assert( y.size() == x.size() );
		
		int n = x.size();
		
		VectorXd dy_dx(n);
		
		dy_dx(0) = ( y(1) - y(0) ) / ( x(1) - x(0) );	// Forward difference.
		
		for ( int i = 1; i < n - 1; ++i ) {
			dy_dx(i) = ( y(i + 1) - y(i - 1) ) / ( x(i + 1) - x(i - 1) );	// Central difference.
		}
		
		dy_dx(n - 1) = ( y(n - 1) - y(n - 2) ) / ( x(n - 1) - x(n - 2) );	// Backward difference.
		
		return dy_dx;
	}
	
	double interp1(const VectorXd & x, const VectorXd & y, const double & xNew)
	{
		assert( x.size() == y.size() );
		assert( x == sort(x) );	// Initial grid must be monotonic.
		
		if ( xNew < x.minCoeff() || xNew > x.maxCoeff() ) {	// New point cannot be external to the initial grid.
			return std::numeric_limits<double>::quiet_NaN();
		}
		
		// Index of the last element in "x" lower than xNew.s: "x(index + 1)" is greater or equal than "xNew".
		int index = std::lower_bound(x.data(), x.data() + x.size(), xNew) - x.data() - 1;
		
		if ( x(index + 1) == xNew ) {	// If "xNew" belongs to the original grid.
			return y(index + 1);
		}
		
		return (double) ( (xNew - x(index + 1)) * y(index) - (xNew - x(index)) * y(index + 1) ) / (x(index) - x(index + 1));
	}
	
	VectorXd interp1(const VectorXd & x, const VectorXd & y, const VectorXd & xNew)
	{
		assert( x.size() == y.size() );
		assert( x == sort(x) );	// Initial grid must be monotonic.
		
		VectorXd yNew = VectorXd::Zero( xNew.size() );
		
		for ( int i = 0; i < yNew.size(); ++i ) {
			yNew(i) = interp1(x, y, xNew(i));
		}
		
		return yNew;
	}
	
	double error_L2(const VectorXd & interp, const VectorXd & simulated,
	                const VectorXd & V, const double & V_shift)
	{
		// Number of not-NaN values.
		unsigned nNotNaN = std::min( (interp   .array() != std::numeric_limits<double>::quiet_NaN()).count(),
		                             (simulated.array() != std::numeric_limits<double>::quiet_NaN()).count()
		                           );
		                           
		VectorXd         V_centered = VectorXd::Zero( nNotNaN );
		VectorXd simulated_centered = VectorXd::Zero( nNotNaN );
		VectorXd    interp_centered = VectorXd::Zero( nNotNaN );
		
		unsigned k = 0;
		
		for ( unsigned i = 0; i < nNotNaN; ++i ) {
			if ( !std::isnan(interp(i)) && !std::isnan(simulated(i)) ) {
				V_centered(k) = V(i);
				simulated_centered(k) = simulated(i);
				interp_centered(k) = interp(i);
				++k;
			}
		}
		
		return trapz(V_centered.array() - V_shift, (interp_centered - simulated_centered).array().square().matrix());
	}
}
