/* C++11 */

/**
 * @file   solvers.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "solvers.h"

PdeSolver1D::PdeSolver1D(VectorXr & mesh)
    : mesh_(mesh), nNodes_(mesh_.size()) {}

Bim1D::Bim1D(VectorXr & mesh)
    : PdeSolver1D(mesh) {}

VectorXr Bim1D::log_mean(const VectorXr & x1, const VectorXr & x2)
{
    assert( x1.size() == x2.size() );
    assert( x1.minCoeff() >= 0.0 );
    assert( x2.minCoeff() >= 0.0 );
    
    VectorXr log_mean = VectorXr::Zero( x1.size() );
    
    for ( Index i = 0; i < log_mean.size(); ++i )
    {
        if ( x1(i) == 0.0 || x2(i) == 0.0 )
        {
            log_mean(i) = 0.0;
        }
        else if ( x1(i) == x2(i) )
        {
            log_mean(i) = x1(i);
        }
        else if ( std::abs( x2(i) - x1(i) ) < 100 * std::numeric_limits<Real>::epsilon() )    // Small values.
        {
            log_mean(i) = 0.5 * ( x1(i) + x2(i) );
        }
        else
        {
            log_mean(i) = ( x2(i) - x1(i) ) / std::log( x2(i) / x1(i) );
        }
    }
    
    return log_mean;
}

std::pair<VectorXr, VectorXr> Bim1D::bernoulli(const VectorXr & x)
{
    VectorXr bp = VectorXr::Zero( x.size() );
    VectorXr bn = VectorXr::Zero( x.size() );
    
    Real lim = 1.0e-2;
    
    VectorXr ax = x.cwiseAbs();
    
    for ( Index i = 0; i < x.size(); ++i )
    {
        if ( x(i) == 0.0 )
        {
            bp(i) = 1.0;
            bn(i) = 1.0;
        }
        else if ( ax(i) > 80.0 )    // Asymptotics.
        {
            if ( x(i) > 0.0 )
            {
                bp(i) = 0.0 ;
                bn(i) = x(i);
            }
            else
            {
                bp(i) = -x(i);
                bn(i) = 0.0  ;
            }
        }
        else if ( ax(i) <= 800.0 && ax(i) > lim )    // Intermediate values.
        {
            bp(i) = x(i) / ( std::exp(x(i)) - 1.0 );
            bn(i) = x(i) + bp(i);
        }
        else    //  if ( ax(i) <= lim ). Small values: Taylor expansion.
        {
            Real  j = 1.0;
            Real fp = 1.0;
            Real fn = 1.0;
            Real df = 1.0;
            int segno = 1  ;
            
            while ( std::abs(df) > 1.0e-16 )
            {
                j += 1.0;
                segno = -segno;
                df *= x(i) / j;
                fp += df;
                fn += segno * df;
            }
            
            bp(i) = 1.0 / fp;
            bn(i) = 1.0 / fn;
        }
    }
    
    return std::make_pair(bp, bn);
}

void Bim1D::assembleAdvDiff(const VectorXr & alpha, const VectorXr & gamma, const VectorXr & eta, const VectorXr & beta)
{
    assert( alpha.size() == nNodes_ - 1 );
    assert( gamma.size() == nNodes_     );
    assert( eta.  size() == nNodes_     );
    
    assert( beta.size() == 1 || beta.size() == nNodes_ || beta.size() == nNodes_ - 1 );
    
    VectorXr area_k = mesh_.segment(1, mesh_.size() - 1) - mesh_.segment(0, mesh_.size() - 1);
    
    VectorXr v_k = VectorXr::Zero( nNodes_ - 1 );
    
    if ( beta.size() == 1)
    {
        v_k.fill(0.0);
    }
    else if ( beta.size() == nNodes_ - 1 )
    {
        v_k = beta.cwiseProduct(area_k);
    }
    else    // if ( beta.size() == nNodes_ )
    {
        v_k = beta.segment(1, beta.size() - 1) - beta.segment(0, beta.size() - 1);
    }
    
    VectorXr gammaEta_k = VectorXr::Zero( nNodes_ - 1 );
    {
        VectorXr temp = gamma.cwiseProduct(eta);
        gammaEta_k = log_mean(temp.segment(0, temp.size() - 1), temp.segment(1, temp.size() - 1));
    }
    
    VectorXr eta_k = log_mean(eta.segment(0, eta.size() - 1), eta.segment(1, eta.size() - 1));
    VectorXr  dEta = eta.segment(1, eta.size() - 1) - eta.segment(0, eta.size() - 1);
    VectorXr   c_k = alpha.cwiseProduct(gammaEta_k).cwiseProduct(eta_k).cwiseQuotient(area_k);
    
    std::pair<VectorXr, VectorXr> bp_bn_k = bernoulli( (v_k - dEta).cwiseQuotient(eta_k) );
    VectorXr & bp = bp_bn_k.first ;
    VectorXr & bn = bp_bn_k.second;
    
    // Assembly
    AdvDiff_.resize( nNodes_, nNodes_ );
    
    for ( Index i = 0; i < AdvDiff_.rows() - 1; ++i )
    {
        AdvDiff_.insert(i + 1,  i ) = - c_k(i) * bp(i);    // Sub-diagonal.
        AdvDiff_.insert( i , i + 1) = - c_k(i) * bn(i);    // Super-diagonal.
    }
    
    AdvDiff_.insert(0, 0) = c_k(0) * bn(0);
    
    for ( Index i = 1; i < AdvDiff_.rows() - 1; ++i )
    {
        AdvDiff_.insert(i, i) = c_k(i) * bn(i) + c_k(i - 1) * bn(i - 1);    // Main diagonal.
    }
    
    AdvDiff_.insert(AdvDiff_.rows() - 1, AdvDiff_.cols() - 1) = c_k(c_k.size() - 1) * bp(bp.size() - 1);
    
    return;
}

void Bim1D::assembleStiff(const VectorXr & eps, const VectorXr & kappa)
{
    assembleAdvDiff( eps, kappa, VectorXr::Ones( nNodes_ ), VectorXr::Zero( 1 ) );
    
    Stiff_ = AdvDiff_;
    
    return;
}

void Bim1D::assembleMass(const VectorXr & delta, const VectorXr & zeta)
{
    assert( delta.size() == nNodes_ - 1 );
    assert( zeta .size() == nNodes_     );
    
    VectorXr h = delta.cwiseProduct( mesh_.segment(1, mesh_.size() - 1) - mesh_.segment(0, mesh_.size() - 1) );
    
    // Assembly
    Mass_.resize( nNodes_, nNodes_ );
    
    Mass_.insert(0, 0) = zeta(0) * 0.5 * h(0);
    
    for ( Index i = 1; i < Mass_.rows() - 1; ++i )
    {
        Mass_.insert(i, i) = zeta(i) * 0.5 * ( h(i - 1) + h(i) );
    }
    
    Mass_.insert(Mass_.rows() - 1, Mass_.cols() - 1) = zeta(zeta.size() - 1) * 0.5 * h(h.size() - 1);
    
    return;
}

NonLinearPoisson1D::NonLinearPoisson1D(const PdeSolver1D & solver, const Index & maxIterationsNo, const Real & tolerance)
    : solver_(solver), maxIterationsNo_(maxIterationsNo), tolerance_(tolerance), qTot_(0.0), cTot_(0.0)/*, cTot_n_(0.0) */
{
    assert( maxIterationsNo_ > 0   );
    assert( tolerance_       > 0.0 );
}

void NonLinearPoisson1D::apply(const VectorXr & mesh, const VectorXr & init_guess, Charge & charge_fun)
{
    assert( init_guess    .size() == mesh.size() );
    assert( solver_.Stiff_.rows() == mesh.size() );
    assert( solver_.Stiff_.cols() == mesh.size() );
    assert( solver_.Mass_ .rows() == mesh.size() );
    assert( solver_.Mass_ .cols() == mesh.size() );
    
    phi_  = VectorXr::Zero( mesh.size() );
    norm_ = VectorXr::Zero( maxIterationsNo_ );
    
    phi_ = init_guess;
    VectorXr phiOld = phi_;
    
    VectorXr  charge = VectorXr::Zero( mesh.size() );
    VectorXr dcharge = VectorXr::Zero( mesh.size() );
    
    SparseXr Jac(solver_.Stiff_.rows(), solver_.Stiff_.rows());
    
    SimplicialLDLT<SparseXr> systemSolver;    // Initialize system solver.
    
    // Newton loop.
    {
        Index k = 0;
        
        for ( ; k < maxIterationsNo_; ++k )
        {
            phiOld.segment(1, phiOld.size() - 2) = phi_.segment(1, phiOld.size() - 2);
            
            charge  = charge_fun. charge(phiOld);
            dcharge = charge_fun.dcharge(phiOld);
            
            // System assembly.
            Jac = computeJac(dcharge);
            
            systemSolver.compute( (SparseXr) Jac.block(1, 1, Jac.rows() - 2, Jac.cols() - 2) );
            
            VectorXr res = (solver_.Stiff_ * phiOld - solver_.Mass_ * charge).segment(1, phiOld.size() - 2);
            
            VectorXr dphi = - systemSolver.solve(res);
            
            /* for ( Index i = 0; i < dphi.size(); ++i )      // Tolerance cut-off.
            {
                if ( dphi(i) > tolerance_ )
                {
                    dphi(i) = tolerance_;
                }
                else if ( dphi(i) < - tolerance_ )
                {
                    dphi(i) = - tolerance_;
                }
            } */
            
            // Newton step.
            phi_.segment(1, phi_.size() - 2) += dphi;    // Dirichlet conditions on boundary.
            
            norm_(k) = dphi.cwiseAbs().maxCoeff();
            
            if ( norm_(k) < tolerance_ )
            {
                break;
            }
        }
        
        if ( k != maxIterationsNo_ - 1 )
        {
            norm_.conservativeResize(k + 1);
        }
    }
    
    for ( Index i = 0; i < phiOld.size(); ++i )
    {
        if ( solver_.Stiff_.coeff(solver_.Stiff_.rows() - 1, i) != 0.0 )
        {
            qTot_ += solver_.Stiff_.coeff(solver_.Stiff_.rows() - 1, i) * phiOld(i);
        }
    }
    
    qTot_ -= solver_.Mass_.coeff(solver_.Mass_.rows() - 1, solver_.Mass_.cols() - 1) * charge(charge.size() - 1);
    
    dcharge = charge_fun.dcharge(phi_);
    
    // System assembly.
    Jac    = computeJac(dcharge);
    
    systemSolver.compute( (SparseXr) Jac.block(1, 1, Jac.rows() - 2, Jac.cols() - 2) );
    
    VectorXr u = VectorXr::LinSpaced(phi_.size(), 0, 1);
    
    // Constant term: b = - Jac(2:end-1, [1 end]) * u([1 end]');
    VectorXr b = VectorXr::Zero( phi_.size() - 2 );
    
    for ( Index i = 1; i < phi_.size() - 1; ++i )
    {
        if ( Jac.coeff(i, 0) != 0.0 || Jac.coeff(i, Jac.cols() - 1) != 0.0 )
        {
            b(i - 1) = Jac.coeff(i, 0) * u(0) + Jac.coeff(i, Jac.cols() - 1) * u(u.size() - 1);
        }
    }
    
    u.segment(1, u.size() - 2) = - systemSolver.solve(b);
    
    for ( Index i = 0; i < phi_.size(); ++i )
    {
        if ( Jac.coeff(Jac.rows() - 1, i) != 0.0 )
        {
            cTot_   += Jac.coeff(Jac.rows() - 1, i) * u(i);
        }
        
        /* if ( Jac.coeff(0, i) != 0.0 )
             {
                cTot_n_ += Jac.coeff(0, i) * u(i);
             } */
    }
    
    // cTot_n_ += cTot_;
}

SparseXr NonLinearPoisson1D::computeJac(const VectorXr & x) const
{
    assert( x.size() == solver_.Stiff_.rows() );
    
    SparseXr Jac( x.size(), x.size() );
    
    Jac = solver_.Stiff_;
    
    for ( Index i = 0; i < Jac.rows(); ++i )
    {
        if ( solver_.Mass_.coeff(i, i) != 0.0 )
        {
            Jac.coeffRef(i, i) -= solver_.Mass_.coeff(i, i) * x(i);
        }
    }
    
    return Jac;
}
