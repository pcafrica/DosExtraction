/* C++11 */

/**
 * @file   charge.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "charge.h"

using namespace constants;

Charge::Charge (const ParamList & params, const QuadratureRule & rule)
  : params_ (params), rule_ (rule) {}

GaussianCharge::GaussianCharge (const ParamList & params,
                                const QuadratureRule & rule)
  : Charge (params, rule) {}

Real
GaussianCharge::n_approx (const Real & phi, const Real & N0,
                               const Real & sigma) const
{
  Real n = 0.0;

  #pragma omp parallel for default(shared) reduction(+: n)
  for (Index i = 0; i < rule_.nNodes_; ++i)
    n += rule_.weights_ (i) * N0 / SQRT_PI /
      (1.0 + std::exp ((SQRT_2 * sigma * rule_.nodes_ (i) - Q * phi) /
                       (K_B * params_.T_)));

  return n;
}

Real
GaussianCharge::dn_approx (const Real & phi, const Real & N0,
                           const Real & sigma) const
{
  Real dn = 0.0;

  #pragma omp parallel for default(shared) reduction(+: dn)
  for (Index i = 0; i < rule_.nNodes_; ++i)
    dn += - Q * rule_.weights_ (i) * N0 * SQRT_2 /
      (sigma * SQRT_PI) * rule_.nodes_ (i) /
      (1.0 + std::exp ((SQRT_2 * sigma * rule_.nodes_
                        (i) - Q * phi) / (K_B * params_.T_)));
  
  return dn;
}

VectorXr
GaussianCharge::charge (const VectorXr & phi) const
{
  VectorXr charge = VectorXr::Zero (phi.size());

  for (Index i = 0; i < phi.size(); ++i)
    {
      charge (i) =
        - Q * n_approx (phi (i), params_.N0_,
                        params_.sigma_);

      if (params_.N0_2_ > 0.0)
        charge (i) +=
          - Q * n_approx (phi (i) + params_.shift_2_,
                          params_.N0_2_, params_.sigma_2_);


      if (params_.N0_3_ > 0.0)
        charge (i) +=
          - Q * n_approx (phi (i) + params_.shift_3_,
                          params_.N0_3_, params_.sigma_3_);

      if (params_.N0_4_ > 0.0)
        charge (i) +=
          - Q * n_approx (phi (i) + params_.shift_4_,
                          params_.N0_4_, params_.sigma_4_);

    }

  return charge;
}

VectorXr
GaussianCharge::dcharge (const VectorXr & phi) const
{
  VectorXr dcharge = VectorXr::Zero (phi.size());

  for (Index i = 0; i < phi.size(); ++i)
    {
      dcharge (i) = - Q * dn_approx (phi (i), params_.N0_,
                                     params_.sigma_);

      if (params_.N0_2_ > 0.0)
        dcharge (i) +=
          - Q * dn_approx (phi (i) + params_.shift_2_,
                           params_.N0_2_, params_.sigma_2_);

      if (params_.N0_3_ > 0.0)
        dcharge (i) +=
          - Q * dn_approx (phi (i) + params_.shift_3_,
                           params_.N0_3_, params_.sigma_3_);

      if (params_.N0_4_ > 0.0)
        dcharge (i) +=
          - Q * dn_approx (phi (i) + params_.shift_4_,
                           params_.N0_4_, params_.sigma_4_);

      if (dcharge (i) > - std::exp (-20.0))
        dcharge (i) = - std::exp (-20.0);
    }

  return dcharge;
}

ExponentialCharge::ExponentialCharge (const ParamList & params,
                                      const QuadratureRule & rule)
  : Charge (params, rule) {}

Real
ExponentialCharge::n_approx (const Real & phi, const Real & N0,
                             const Real & lambda) const
{
  Real n = 0.0;

  #pragma omp parallel for default(shared) reduction(+: n)
  for (Index i = 0; i < rule_.nNodes_; ++i)
    n += rule_.weights_ (i) * N0 /
      (1.0 + std::exp ((lambda * rule_.nodes_ (i) - Q * phi) /
                       (K_B * params_.T_)));

  return n;
}

Real
ExponentialCharge::dn_approx (const Real & phi, const Real & N0,
                              const Real & lambda) const
{
  Real dn = 0.0;

  #pragma omp parallel for default(shared) reduction(+: dn)
  for (Index i = 0; i < rule_.nNodes_; ++i)
    dn +=
      - Q * rule_.weights_ (i) * N0 / (lambda) * rule_.nodes_ (i) /
      (1.0 + std::exp ((lambda * rule_.nodes_ (i) - Q * phi) /
                       (K_B * params_.T_)));

  return dn;
}

VectorXr
ExponentialCharge::charge (const VectorXr & phi) const
{
  VectorXr charge = VectorXr::Zero (phi.size());

  for (Index i = 0; i < phi.size(); ++i)
    charge (i) =
      - Q * n_approx (phi (i), params_.N0_exp_,
                      params_.lambda_exp_);
  return charge;
}

VectorXr
ExponentialCharge::dcharge (const VectorXr & phi) const
{
  VectorXr dcharge = VectorXr::Zero (phi.size());

  for (Index i = 0; i < phi.size(); ++i)
    dcharge (i) =
      - Q * dn_approx (phi (i), params_.N0_exp_,
                       params_.lambda_exp_);

  return dcharge;
}

