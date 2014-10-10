/* C++11 */

/**
 * @file   quadratureRule.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "quadratureRule.h"

using namespace constants;

QuadratureRule::QuadratureRule(const Index & nNodes)
  : nNodes_(nNodes)
{
  assert( nNodes_ >= 1 );
  
  nodes_  .resize( nNodes_ );
  weights_.resize( nNodes_ );
}

GaussHermiteRule::GaussHermiteRule(const Index & nNodes)
  : QuadratureRule(nNodes) {}

void GaussHermiteRule::apply()
{
  apply_iterative_algorithm();    // Using default parameters for maximum iterations number and tolerance.
}

void GaussHermiteRule::apply(const GetPot & config)
{
  apply_iterative_algorithm( config("QuadratureRule/maxIterationsNo", 1000), config("QuadratureRule/tolerance", 1.0e-14) );
}

void GaussHermiteRule::apply_iterative_algorithm(const Index & maxIterationsNo, const Real & tolerance)
{
  assert( maxIterationsNo > 0   );
  assert( tolerance       > 0.0 );
  
  Real p1 = 0.0, p2 = 0.0, temp = 0.0, dp = 0.0;
  Real z = 0.0, zOld = 0.0;
  
  // The roots are symmetric about the origin: finding only half of them is needed.
  for ( Index i = 0; i < (nNodes_ + 1) / 2; ++i )
    {
      // Initial guesses for the largest roots.
      if ( i == 0 )
        {
          z = std::sqrt(2.0 * nNodes_ + 1.0) - 1.85575 * std::pow(2.0 * nNodes_ + 1.0, - 0.16667);
        }
      else if ( i == 1 )
        {
          z -= 1.14 * std::pow(nNodes_, 0.426) / z;
        }
      else if ( i == 2 )
        {
          z = 1.86 * z - 0.86 * nodes_(nNodes_ - 1);
        }
      else if ( i == 3 )
        {
          z = 1.91 * z - 0.91 * nodes_(nNodes_ - 2);
        }
      else
        {
          z = 2.0 * z - nodes_(nNodes_ - i + 1);
        }
        
      Index j = 0;
      
      for ( ; j < maxIterationsNo; ++j )    // Refinement by Newton's method.
        {
          p1 = PI_M4;
          p2 = 0.0;
          
          for ( Index k = 0; k < nNodes_; ++k )    // Loop up the recurrence relation to evaluate the Hermite polynomial at "z".
            {
              temp = p2;
              p2 = p1;
              p1 = z * std::sqrt( 2.0 / (k + 1.0) ) * p2 - std::sqrt( k / (k + 1.0) ) * temp;
            }
            
          // "p1" is now the desidered Hermite polynomial. Computing its derivative "dp" is needed.
          dp = std::sqrt(2.0 * nNodes_) * p2;
          
          // Newton step
          zOld = z;
          z   -= p1 / dp;
          
          if ( std::abs(z - zOld) <= tolerance )
            {
              break;
            }
        }
        
      if ( j == maxIterationsNo )
        {
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
  if ( nNodes_ == 1 )
    {
      nodes_  .fill(0.0);
      weights_.fill(PI) ;
      return;
    }
    
  MatrixXr Jac = MatrixXr::Zero( nNodes_, nNodes_ );    // Jacobi matrix.
  VectorXr k = VectorXr::LinSpaced( nNodes_, 1, nNodes_ );
  VectorXr v = std::sqrt(0.5) * k.cwiseSqrt();
  
  for ( Index i = 0; i < Jac.rows() - 1; ++i )
    {
      Jac(i + 1,   i  ) = v(i);    // Sub-diagonal.
      Jac(  i  , i + 1) = v(i);    // Super-diagonal.
    }
    
  SelfAdjointEigenSolver<MatrixXr> eigSolver(Jac);
  
  nodes_ = eigSolver.eigenvalues();
  
  MatrixXr EigVec = eigSolver.eigenvectors();
  
  VectorXr norm2 = (EigVec.transpose() * EigVec ).diagonal().cwiseSqrt();    // Normalize weights.
  
  weights_ = ( SQRT_PI * EigVec.row(0).cwiseProduct(EigVec.row(0)).transpose() ).cwiseQuotient(norm2);    // SQRT_PI = beta0;
  
  return;
}

GaussLaguerreRule::GaussLaguerreRule(const Index & nNodes)
  : QuadratureRule(nNodes) {}

Real GaussLaguerreRule::log_gamma(const Real & x)
{
  assert( x > 0.0 );
  
  static const Real coeff[14] = { 57.1562356658629235, -59.5979603554754912,
                                  14.1360979747417471, -0.491913816097620199,
                                  0.339946499848118887e-4, 0.465236289270485756e-4,
                                  -0.983744753048795646e-4, 0.158088703224912494e-3,
                                  -0.210264441724104883e-3, 0.217439618115212643e-3,
                                  -0.164318106536763890e-3, 0.844182239838527433e-4,
                                  -0.261908384015814087e-4, 0.368991826595316234e-5
                                };
                                
  Real temp = x + 5.24218750000000000;    // x + 671/128
  
  temp = (x + 0.5) * std::log(temp) - temp;
  
  Real aux = 0.999999999999997092;
  
  for ( Index j = 0; j < 14; ++j )
    {
      aux += coeff[j] / (x + j);
    }
    
  return temp + std::log(2.5066282746310005 * aux / x);
}

void GaussLaguerreRule::apply()
{
  apply_iterative_algorithm();    // Using default parameters for maximum iterations number and tolerance.
}

void GaussLaguerreRule::apply(const GetPot & config)
{
  apply_iterative_algorithm( config("QuadratureRule/maxIterationsNo", 1000), config("QuadratureRule/tolerance", 1.0e-14) );
}

void GaussLaguerreRule::apply_iterative_algorithm(const Index & maxIterationsNo, const Real & tolerance)
{
  assert( maxIterationsNo > 0   );
  assert( tolerance       > 0.0 );
  
  Real p1 = 0.0, p2 = 0.0, temp = 0.0, dp = 0.0;
  Real z = 0.0, zOld = 0.0;
  
  // Loop over the desired roots.
  for ( Index i = 0; i < nNodes_; ++i )
    {
      // Initial guesses for the largest roots.
      if ( i == 0 )
        {
          z = 3.0 / ( 1.0 + 2.4 * nNodes_ );
        }
      else if ( i == 1 )
        {
          z += 15.0 / ( 1.0 + 2.5 * nNodes_ );
        }
      else
        {
          z += ( (1.0 + 2.55 * (i - 1)) / (1.9 * (i - 1)) ) *  (z - nodes_(i - 2));
        }
        
      Index j = 0;
      
      for ( ; j < maxIterationsNo; ++j )    // Refinement by Newton's method.
        {
          p1 = 1.0;
          p2 = 0.0;
          
          for ( Index k = 0; k < nNodes_; ++k )    // Loop up the recurrence relation to evaluate the Hermite polynomial at "z".
            {
              temp = p2;
              p2 = p1;
              p1 = ( (2.0 * k + 1.0 - z) * p2 - k * temp ) / ( k + 1.0 );
            }
            
          // "p1" is now the desidered Hermite polynomial. Computing its derivative "dp" is needed.
          dp = nNodes_ * (p1 - p2) / z;
          
          // Newton step
          zOld = z;
          z   -= p1 / dp;
          
          if ( std::abs(z - zOld) <= tolerance )
            {
              break;
            }
        }
        
      if ( j == maxIterationsNo )
        {
          throw std::runtime_error("ERROR: GaussLaguerreRule::apply_iterative_algorithm() didn't reach convergence.");
        }
        
      nodes_(i) = z;
      
      weights_(i) = - 1 / (dp * nNodes_ * p2);
      
    }
    
  return;
}

void GaussLaguerreRule::apply_using_eigendecomposition()
{
  if ( nNodes_ == 1 )
    {
      nodes_  .fill(0.0);
      weights_.fill(1.0);
      return;
    }
    
  MatrixXr Jac = MatrixXr::Zero( nNodes_, nNodes_ );    // Jacobi matrix.
  VectorXr k = VectorXr::LinSpaced( nNodes_, 1, nNodes_ );
  
  for ( Index i = 0; i < Jac.rows() - 1; ++i )
    {
      Jac(i + 1,   i  ) = k(i);    // Sub-diagonal.
      Jac(  i  , i + 1) = k(i);    // Super-diagonal.
    }
    
  for ( Index i = 0; i < Jac.rows(); ++i )
    {
      Jac(i, i) = 2 * k(i) - 1;    // Main diagonal.
    }
    
  SelfAdjointEigenSolver<MatrixXr> eigSolver(Jac);
  
  nodes_ = eigSolver.eigenvalues();
  
  MatrixXr EigVec = eigSolver.eigenvectors();
  
  VectorXr norm2 = (EigVec.transpose() * EigVec ).diagonal().cwiseSqrt();    // Normalize weights.
  
  weights_ = ( 1 * EigVec.row(0).cwiseProduct(EigVec.row(0)).transpose() ).cwiseQuotient(norm2);    // 1 = beta0;
  
  return;
}
