/* C++11 */

/**
 * @file   quadratureRule.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 */

#include "quadratureRule.h"

using namespace constants;

QuadratureRule::QuadratureRule(const Index & nNodes)
	: nNodes_(nNodes)
{
	assert(nNodes_ >= 1);
	
	nodes_  .resize( nNodes_ );
	weights_.resize( nNodes_ );
}

GaussHermiteRule::GaussHermiteRule(const Index & nNodes)
	: QuadratureRule(nNodes) {}

void GaussHermiteRule::apply()
{
	apply_iterative_algorithm();	// Using default parameters for maximum iterations number and tolerance.
}

void GaussHermiteRule::apply(const GetPot & config)
{
	apply_iterative_algorithm( config("GaussHermite/maxIterationsNo", 1000), config("GaussHermite/tolerance", 1.0e-14) );
}

void GaussHermiteRule::apply_iterative_algorithm(const Index & maxIterationsNo, const Real & tolerance)
{
	Real p1, p2, temp, dp;
	Real z, zOld;
	
	for ( Index i = 0; i < (nNodes_ + 1) / 2; ++i ) {	// The roots are symmetric about the origin,
		// so only finding half of them is needed.
		// Good guesses at initial values for the most largest roots.
		if ( i == 0 ) {
			z = std::sqrt(2.0 * nNodes_ + 1.0) - 1.85575 * std::pow(2.0 * nNodes_ + 1.0, - 0.16667);
		} else if ( i == 1 ) {
			z -= 1.14 * std::pow(nNodes_, 0.426) / z;
		} else if ( i == 2 ) {
			z = 1.86 * z - 0.86 * nodes_(nNodes_ - 1);
		} else if ( i == 3 ) {
			z = 1.91 * z - 0.91 * nodes_(nNodes_ - 2);
		} else {
			z = 2.0 * z - nodes_(nNodes_ - i + 1);
		}
		
		Index j = 0;
		
		for ( ; j < maxIterationsNo; ++j ) {	// Refinement by Newton's method.
			p1 = PI_M4;
			p2 = 0.0;
			
			for ( Index k = 0; k < nNodes_; ++k ) {	// Loop up the recurrence relation to evaluate the Hermite polynomial at "z".
				temp = p2;
				p2 = p1;
				p1 = z * std::sqrt( 2.0 / (k + 1.0) ) * p2 - std::sqrt( k / (k + 1.0) ) * temp;
			}
			
			// "p1" is now the desidered Hermite polynomial. Computing its derivative "dp" is needed.
			dp = std::sqrt(2.0 * nNodes_) * p2;
			
			// Newton step
			zOld = z;
			z  = zOld - p1 / dp;
			
			if ( std::abs(z - zOld) <= tolerance ) {
				break;
			}
		}
		
		if ( j >= maxIterationsNo ) {
			throw std::runtime_error("ERROR: GaussHermiteRule::apply_iterative_algorithm() didn't reach convergence.");
		}
		
		nodes_(i)               = -z;
		nodes_(nNodes_ - i - 1) =  z;
		
		weights_(i)               = 2.0 / (dp * dp);
		weights_(nNodes_ - i - 1) = weights_(i)   ;
	}
	
	return;
}

void GaussHermiteRule::apply_using_eigendecomposition()
{
	if (nNodes_ == 1) {
		nodes_  .fill(0.0);
		weights_.fill(PI) ;
		return;
	}
	
	MatrixXr Jac = MatrixXr::Zero( nNodes_, nNodes_ );	// Jacobi matrix.
	VectorXr k = VectorXr::LinSpaced( nNodes_, 1, nNodes_);
	VectorXr v = std::sqrt(0.5) * k.cwiseSqrt();
	
	for ( Index i = 0; i < Jac.rows() - 1; ++i ) {
		Jac(i + 1,  i ) = v(i);	// Sub-diagonal.
		Jac( i , i + 1) = v(i);	// Super-diagonal.
	}
	
	SelfAdjointEigenSolver<MatrixXr> eigSolver(Jac);
	
	nodes_ = eigSolver.eigenvalues();
	
	MatrixXr EigVec = eigSolver.eigenvectors();
	
	VectorXr norm2 = (EigVec.transpose() * EigVec ).diagonal().cwiseSqrt();	// To normalize weights.
	
	weights_ = ( SQRT_PI * EigVec.row(0).cwiseProduct(EigVec.row(0)).transpose() ).cwiseQuotient(norm2);	// sqrt(pi) = beta0;
	
	return;
}
