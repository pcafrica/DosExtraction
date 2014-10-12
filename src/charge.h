/* C++ */

/**
 * @file   charge.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Classes for computing total electric charge density.
 *
 */

#ifndef CHARGE_H
#define CHARGE_H

#include "paramList.h"
#include "quadratureRule.h"
#include "typedefs.h"

/**
 * @class Charge
 *
 * @brief Abstract class providing methods to calculate total electric charge density (the rhs in the Poisson equation).
 *
 */
class Charge
{
    public:
        /**
         * @brief Default constructor (deleted since it is required to specify a @ref ParamList and a @ref QuadratureRule).
         */
        Charge() = delete;
        /**
         * @brief Constructor.
         * @param[in] params : a list of simulation parameters;
         * @param[in] rule   : a quadrature rule.
         */
        Charge(const ParamList &, const QuadratureRule &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~Charge() = default;
        
        /**
         * @brief Compute the total charge density.
         * @param[in] phi : the electric potential @f$ \varphi @f$.
         * @returns the total charge density @f$ q(\varphi) \left[ C \cdot m^{-3} \right] @f$.
         */
        virtual VectorXr  charge(const VectorXr & phi) = 0;
        /**
         * @brief Compute the derivative of the total charge density with respect to the electric potential.
         * @param[in] phi : the electric potential @f$ \varphi @f$.
         * @returns the derivative: @f$ \frac{\mathrm{d}q(\varphi)}{\mathrm{d}\varphi} \left[ C \cdot m^{-3} \cdot V^{-1} \right] @f$.
         */
        virtual VectorXr dcharge(const VectorXr & phi) = 0;
        
    protected:
        const ParamList      & params_;    /**< @brief Parameter list handler. */
        const QuadratureRule & rule_  ;    /**< @brief Quadrature rule handler. */
};

/**
 * @class GaussianCharge
 *
 * Provide methods to compute total electric charge and its derivative under the hypothesis that Density of States
 * is a linear combination of multiple gaussians, whose parameters are read from a @ref ParamList object, of the form:
 * @f[ \frac{N_0}{\sqrt{2\pi\sigma^2}}\exp\left(-\frac{\left(\cdot\right)^2}{2\sigma^2}\right) ~ . @f]
 *
 * @brief Class derived from @ref Charge, under the hypothesis that Density of States is a combination of gaussians.
 *
 */
class GaussianCharge : public Charge
{
    public:
        /**
         * @brief Default constructor (deleted since it is required to specify a @ref ParamList and a @ref QuadratureRule).
         */
        GaussianCharge() = delete;
        /**
         * @brief Constructor.
         * @param[in] params : a list of simulation parameters;
         * @param[in] rule   : a quadrature rule.
         */
        GaussianCharge(const ParamList &, const QuadratureRule &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~GaussianCharge() = default;
        
        virtual VectorXr  charge(const VectorXr &) override;
        virtual VectorXr dcharge(const VectorXr &) override;
        
    private:
        /**
         * @brief Compute electrons density (per unit volume).
         * @param[in] phi   : the electric potential @f$ \varphi @f$;
         * @param[in] N0    : the gaussian mean @f$ N_0 @f$;
         * @param[in] sigma : the gaussian standard deviation @f$ \sigma @f$.
         * @returns the electrons density @f$ n(\varphi) \left[ m^{-3} \right] @f$.
         */
        Real  n_approx(const Real &, const Real &, const Real &) const;
        /**
         * @brief Compute the approximate derivative of electrons density (per unit volume)
         * with respect to the electric potential.
         * @param[in] phi   : the electric potential @f$ \varphi @f$;
         * @param[in] N0    : the gaussian mean @f$ N_0 @f$;
         * @param[in] sigma : the gaussian standard deviation @f$ \sigma @f$.
         * @returns the derivative: @f$ \frac{\mathrm{d}n(\varphi)}{\mathrm{d}\varphi} \left[ m^{-3} \cdot V^{-1} \right] @f$.
         */
        Real dn_approx(const Real &, const Real &, const Real &) const;
};

/**
 * @class ExponentialCharge
 *
 * Provide methods to compute total electric charge and its derivative under the hypothesis that Density of States
 * is a single exponential, whose parameter is got by the constructor, of the form:
 * @f[ \frac{N_0}{\lambda}\exp\left(-\frac{\left(\cdot\right)}{\lambda}\right) ~ . @f]
 *
 * @brief Class derived from @ref Charge, under the hypothesis that Density of States is a single exponential.
 *
 */
class ExponentialCharge : public Charge
{
    public:
        /**
         * @brief Default constructor (deleted since it is required to specify a @ref ParamList and a @ref QuadratureRule).
         */
        ExponentialCharge() = delete;
        /**
         * @brief Constructor.
         * @param[in] params : a list of simulation parameters;
         * @param[in] rule : a quadrature rule.
         */
        ExponentialCharge(const ParamList &, const QuadratureRule &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~ExponentialCharge() = default;
        
        virtual VectorXr  charge(const VectorXr &) override;
        virtual VectorXr dcharge(const VectorXr &) override;
        
    private:
        /**
         * @brief Compute electrons density (per unit volume).
         * @param[in] phi    : the electric potential @f$ \varphi @f$;
         * @param[in] N0     : the exponential @f$ N_0 @f$;
         * @param[in] lambda : the exponential @f$ \lambda @f$.
         * @returns the electrons density @f$ n(\varphi) \left[ m^{-3} \right] @f$.
         */
        Real  n_approx(const Real &, const Real &, const Real &) const;
        /**
         * @brief Compute the approximate derivative of electrons density (per unit volume)
         * with respect to the electric potential.
         * @param[in] phi : the electric potential @f$ \varphi @f$;
         * @param[in] N0     : the exponential @f$ N_0 @f$;
         * @param[in] lambda : the exponential @f$ \lambda @f$.
         * @returns the derivative: @f$ \frac{\mathrm{d}n(\varphi)}{\mathrm{d}\varphi} \left[ m^{-3} \cdot V^{-1} \right] @f$.
         */
        Real dn_approx(const Real &, const Real &, const Real &) const;
};

#endif /* CHARGE_H */
