#include "charge.h"

using namespace constants;

Charge::Charge(const ParamList & params, const QuadratureRule & rule)
	: params_(params), rule_(rule) {}

GaussianCharge::GaussianCharge(const ParamList & params, const QuadratureRule & rule)
	: Charge(params, rule) {}

double GaussianCharge::n_approx(const double & phi, const double & N0, const double & sigma) const
{
	double n = 0.0;
	
	#pragma omp parallel for default(shared) reduction(+: n)
	
	for (unsigned i = 0; i < rule_.nNodes_; ++i) {
		n += rule_.weights_(i) * N0 / SQRT_PI / ( 1.0 + std::exp( (SQRT_2 * sigma * rule_.nodes_(i) - Q * phi) / (K_B * T) ) );
	}
	
	return n;
}

double GaussianCharge::dn_approx(const double & phi, const double & N0, const double & sigma) const
{
	double dn = 0.0;
	
	#pragma omp parallel for default(shared) reduction(+: dn)
	
	for ( unsigned i = 0; i < rule_.nNodes_; ++i) {
		dn += rule_.weights_(i) * N0 * SQRT_2 / (sigma * SQRT_PI) * rule_.nodes_(i) / ( 1.0 + std::exp( (SQRT_2 * sigma * rule_.nodes_(i) - Q * phi) / (K_B * T) ) );
	}
	
	return dn;
}

VectorXd GaussianCharge::charge(const VectorXd & phi)
{
	VectorXd charge = VectorXd::Zero( phi.size() );
	
	for ( unsigned i = 0; i < phi.size(); ++i ) {
		charge(i) = - Q * n_approx(phi(i), params_.N0_, params_.sigma_);
		
		if ( params_.N0_2_ > 0.0 ) {
			charge(i) += - Q * n_approx(phi(i) + params_.shift_2_, params_.N0_2_, params_.sigma_2_);
		}
		
		if ( params_.N0_3_ > 0.0 ) {
			charge(i) += - Q * n_approx(phi(i) + params_.shift_3_, params_.N0_3_, params_.sigma_3_);
		}
		
		if ( params_.N0_4_ > 0.0 ) {
			charge(i) += - Q * n_approx(phi(i) + params_.shift_4_, params_.N0_4_, params_.sigma_4_);
		}
	}
	
	return charge;
}

VectorXd GaussianCharge::dcharge(const VectorXd & phi)
{
	VectorXd dcharge = VectorXd::Zero( phi.size() );
	
	for ( unsigned i = 0; i < phi.size(); ++i ) {
		dcharge(i) = Q2 * dn_approx(phi(i), params_.N0_, params_.sigma_);
		
		if ( params_.N0_2_ > 0.0 ) {
			dcharge(i) += Q2 * dn_approx(phi(i) + params_.shift_2_, params_.N0_2_, params_.sigma_2_);
		}
		
		if ( params_.N0_3_ > 0.0 ) {
			dcharge(i) += Q2 * dn_approx(phi(i) + params_.shift_3_, params_.N0_3_, params_.sigma_3_);
		}
		
		if ( params_.N0_4_ > 0.0 ) {
			dcharge(i) += Q2 * dn_approx(phi(i) + params_.shift_4_, params_.N0_4_, params_.sigma_4_);
		}
		
		if ( dcharge(i) > - std::exp(-20) ) {
			dcharge(i) = - std::exp(-20);
		}
	}
	
	return dcharge;
}
