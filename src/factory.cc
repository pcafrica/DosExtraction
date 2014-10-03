/* C++11 */

/**
 * @file   factory.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "factory.h"

Charge * GaussianChargeFactory::BuildCharge(const ParamList & params, const QuadratureRule & rule)
{
  return new GaussianCharge(params, rule);
}

Charge * ExponentialChargeFactory::BuildCharge(const ParamList & params, const QuadratureRule & rule)
{
  return new ExponentialCharge(params, rule);
}

QuadratureRule * GaussHermiteRuleFactory::BuildRule(const Index & nNodes)
{
  return new GaussHermiteRule(nNodes);
}

QuadratureRule * GaussLaguerreRuleFactory::BuildRule(const Index & nNodes)
{
  return new GaussLaguerreRule(nNodes);
}
