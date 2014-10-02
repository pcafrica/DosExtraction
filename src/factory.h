/* C++ */

/**
 * @file   factory.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Abstract factories for the DOS constitutive relation and for the quadrature rule.
 *
 */

#ifndef FACTORY_H
#define FACTORY_H

#include "charge.h"
#include "paramList.h"
#include "quadratureRule.h"

/**
 * @class ChargeFactory
 *
 * @brief Abstract factory to handle at runtime a constitutive relation for the Density of States.
 *
 */
class ChargeFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    ChargeFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~ChargeFactory() = default;
    
    /**
     * @brief Method to build an abstract @ref Charge object.
     * @param[in] params : the list of simulation parameters;
     * @param[in] rule   : a quadrature rule.
     * @returns a pointer to @ref Charge.
     */
    virtual Charge * BuildCharge(const ParamList & params, const QuadratureRule & rule) = 0;
};

/**
 * @class GaussianChargeFactory
 *
 * @brief Concrete factory to handle a multiple gaussians constitutive relation.
 *
 */
class GaussianChargeFactory : public ChargeFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    GaussianChargeFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~GaussianChargeFactory() = default;
    
    /**
     * @brief Method to build a concrete @ref Charge object.
     * @param[in] params : the list of simulation parameters;
     * @param[in] rule   : a quadrature rule.
     * @returns a pointer to @ref GaussianCharge.
     */
    virtual Charge * BuildCharge(const ParamList &, const QuadratureRule &) override;
};

/**
 * @class ExponentialChargeFactory
 *
 * @brief Concrete factory to handle a single exponential constitutive relation.
 *
 */
class ExponentialChargeFactory : public ChargeFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    ExponentialChargeFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~ExponentialChargeFactory() = default;
    
    /**
     * @brief Method to build a concrete @ref Charge object.
     * @param[in] params : the list of simulation parameters;
     * @param[in] rule   : a quadrature rule.
     * @returns a pointer to @ref ExponentialCharge.
     */
    virtual Charge * BuildCharge(const ParamList &, const QuadratureRule &) override;
};

/**
 * @class QuadratureRuleFactory
 *
 * @brief Abstract factory to handle at runtime a quadrature rule.
 *
 */
class QuadratureRuleFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    QuadratureRuleFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~QuadratureRuleFactory() = default;
    
    /**
     * @brief Method to build an abstract @ref QuadratureRule object.
     * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
     * @returns a pointer to @ref QuadratureRule.
     */
    virtual QuadratureRule * BuildRule(const Index & nNodes) = 0;
};

/**
 * @class GaussHermiteRuleFactory
 *
 * @brief Concrete factory to handle a Gauss-Hermite quadrature rule.
 *
 */
class GaussHermiteRuleFactory : public QuadratureRuleFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    GaussHermiteRuleFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~GaussHermiteRuleFactory() = default;
    
    /**
     * @brief Method to build a concrete @ref QuadratureRule object.
     * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
     * @returns a pointer to @ref GaussHermiteRule.
     */
    virtual QuadratureRule * BuildRule(const Index &) override;
};

/**
 * @class GaussLaguerreRuleFactory
 *
 * @brief Concrete factory to handle a Gauss-Laguerre quadrature rule.
 *
 */
class GaussLaguerreRuleFactory : public QuadratureRuleFactory
{
  public:
    /**
     * @brief Default constructor (defaulted).
     */
    GaussLaguerreRuleFactory() = default;
    /**
     * @brief Destructor (defaulted).
     */
    virtual ~GaussLaguerreRuleFactory() = default;
    
    /**
     * @brief Method to build a concrete @ref QuadratureRule object.
     * @param[in] nNodes : the number of nodes to be used for the quadrature rule.
     * @returns a pointer to @ref GaussLaguerreRule.
     */
    virtual QuadratureRule * BuildRule(const Index &) override;
};

#endif /* FACTORY_H */
