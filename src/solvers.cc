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
        else if ( ax(i) <= 80.0 && ax(i) > lim )    // Intermediate values.
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
            
            int sign = 1;
            
            while ( std::abs(df) > 1.0e-16 )
            {
                j += 1.0;
                sign = -sign;
                df *= x(i) / j;
                fp += df;
                fn += sign * df;
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
    
    if ( beta.size() == 1 )
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
        AdvDiff_.insert(i + 1,  i  ) = - c_k(i) * bp(i);    // Sub-diagonal.
        AdvDiff_.insert(  i , i + 1) = - c_k(i) * bn(i);    // Super-diagonal.
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
    assembleAdvDiff( eps, kappa, VectorXr::Ones( nNodes_ ), VectorXr::Zero(1) );
    
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

NonLinearPoisson1D::NonLinearPoisson1D(const ParamList & params, const PdeSolver1D & solver, const Index & maxIterationsNo, const Real & tolerance)
    : params_(params), solver_(solver), maxIterationsNo_(maxIterationsNo), tolerance_(tolerance), PhiBcorr_(0.0), qTot_(0.0), cTot_(0.0)/*, cTot_n_(0.0) */
{
    assert( maxIterationsNo_ > 0   );
    assert( tolerance_       > 0.0 );
}

void NonLinearPoisson1D::apply(const VectorXr & init_guess, const Charge & charge_fun)
{
    assert( init_guess    .size() == solver_.mesh_.size() );
    assert( solver_.Stiff_.rows() == solver_.mesh_.size() );
    assert( solver_.Stiff_.cols() == solver_.mesh_.size() );
    assert( solver_.Mass_ .rows() == solver_.mesh_.size() );
    assert( solver_.Mass_ .cols() == solver_.mesh_.size() );
    
    phi_      = init_guess;
    norm_     = VectorXr::Zero( maxIterationsNo_ );
    
    PhiBcorr_ = 0.0;
    qTot_     = 0.0;
    
    VectorXr phiOld = phi_;
    
    VectorXr  charge = VectorXr::Zero( solver_.mesh_.size() );
    VectorXr dcharge = VectorXr::Zero( solver_.mesh_.size() );
    
    SparseXr Jac(solver_.Stiff_.rows(), solver_.Stiff_.rows());
    
    SparseLU<SparseXr> systemSolver;    // Initialize system solver.
    
    // Newton loop.
    {
        Index k = 0;
        
        for ( ; k < maxIterationsNo_; ++k )
        {
            phiOld = phi_;
            
            charge  = charge_fun. charge(phiOld.array() + constants::V_TH * PhiBcorr_);
            dcharge = charge_fun.dcharge(phiOld.array() + constants::V_TH * PhiBcorr_);
            
            // System assembly.
            VectorXr res = (solver_.Stiff_ * phiOld - solver_.Mass_ * charge);
            
            Real E = res(0) / params_.eps_semic();
            
            const Real coeff = params_.PhiBcoeff();
            const Real f = coeff * coeff * (-E); // Negative field => electron injection.
            
            if (f > 0)
            {
                PhiBcorr_ = std::sqrt(f);
            }
            else
            {
                PhiBcorr_ = f / 4;
            }
            
            Jac = computeJac(dcharge);
            
            systemSolver.compute( (SparseXr) Jac.block(1, 1, Jac.rows() - 2, Jac.cols() - 2) );
            
            VectorXr dphi = - systemSolver.solve(res.segment(1, phiOld.size() - 2));
            
            /*for ( Index i = 0; i < dphi.size(); ++i )    // Damping.
            {
                if ( dphi(i) > tolerance_ )
                {
                    dphi(i) = tolerance_;
                }
                else if ( dphi(i) < - tolerance_ )
                {
                    dphi(i) = - tolerance_;
                }
            }*/
            
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
    
    dcharge = charge_fun.dcharge(phi_.array() + PhiBcorr_);
    
    // Compute total capacitance.
    // System assembly.
    Jac = computeJac(dcharge);
    
    systemSolver.compute( (SparseXr) Jac.block(1, 1, Jac.rows() - 2, Jac.cols() - 2) );
    
    VectorXr u = VectorXr::LinSpaced(phi_.size(), 0, 1);
    
    // Constant term: b = - Jac(2:end-1, [1 end]) * u([1 end]');
    VectorXr b = -(Jac.block(1, 0, Jac.rows() - 2, 1) * u(0) +
                   Jac.block(1, Jac.cols() - 1, Jac.rows() - 2, 1) * u(u.size() - 1));
                   
    u.segment(1, u.size() - 2) = systemSolver.solve(b);
    
    cTot_ = ((VectorXr) Jac.row(Jac.rows() - 1)).dot(u);
    /* cTot_n_ = ((VectorXr) Jac.row(0)).dot(u) + cTot_; */
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
