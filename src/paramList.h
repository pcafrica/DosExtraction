/* C++ */

/**
 * @file   paramList.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Interface to process a list of simulation parameters.
 *
 */

#ifndef PARAMLIST_H
#define PARAMLIST_H

#include "typedefs.h"

/**
 * @class ParamList
 *
 * It can include up to 4 gaussians (later combined to compute total charge) and an exponential.
 *
 * @brief Class providing methods to handle a list of parameters.
 *
 */
class ParamList
{
    public:
        friend class    GaussianCharge;
        friend class ExponentialCharge;
        friend class          DosModel;
        // Now GaussianCharge, ExponentialCharge and DosModel can access private parameters with
        // no need to copy them through getter methods that slow down the program execution.
        
        /**
         * @brief Default constructor (defaulted).
         */
        ParamList() = default;
        /**
         * @brief Explicit conversion constructor.
         * @param[in] list : a row vector containing a list of parameters (for example got by a @ref CsvParser object).
         * Parameters should be sorted in the same order as specified above.
         */
        explicit ParamList(const RowVectorXr &);
        /**
         * @brief Destructor (defaulted).
         */
        virtual ~ParamList() = default;
        
        /**
         * @name Getter methods
         * @{
         */
        inline const Index & simulationNo() const;
        inline const Real  & t_semic()      const;
        inline const Real  & t_ins()        const;
        inline const Real  & eps_semic()    const;
        inline const Real  & eps_ins()      const;
        inline const Real  & Wf()           const;
        inline const Real  & Ea()           const;
        inline const Real  & N0()           const;
        inline const Real  & sigma()        const;
        inline const Real  & N0_2()         const;
        inline const Real  & sigma_2()      const;
        inline const Real  & shift_2()      const;
        inline const Real  & N0_3()         const;
        inline const Real  & sigma_3()      const;
        inline const Real  & shift_3()      const;
        inline const Real  & N0_4()         const;
        inline const Real  & sigma_4()      const;
        inline const Real  & shift_4()      const;
        inline const Real  & N0_exp()       const;
        inline const Real  & lambda_exp()   const;
        inline const Real  & A_semic()      const;
        inline const Real  & C_sb()         const;
        inline const Index & nNodes()       const;
        inline const Index & nSteps()       const;
        inline const Real  & V_min()        const;
        inline const Real  & V_max()        const;
        
        /**
         * @}
         */
        
    private:
        Index simulationNo_;    /**< @brief Simulation number index. */
        Real  t_semic_     ;    /**< @brief Thickness of the semiconductor layer @f$ \left[ m \right] @f$. */
        Real  t_ins_       ;    /**< @brief Thickness of the insulator layer @f$ \left[ m \right] @f$. */
        Real  eps_semic_   ;    /**< @brief Relative electrical permittivity of the semiconductor layer @f$ \left[ ~ \right] @f$. */
        Real  eps_ins_     ;    /**< @brief Relative electrical permittivity of the insulator layer @f$ \left[ ~ \right] @f$. */
        Real  Wf_          ;    /**< @brief Work-function (normalized by @f$ Q @f$) @f$ \left[ V \right] @f$. */
        Real  Ea_          ;    /**< @brief Electron affinity (normalized by @f$ Q @f$) @f$ \left[ V \right] @f$. */
        Real  N0_          ;    /**< @brief 1st gaussian @f$ N_0 \left[ m^{-3} \right] @f$. */
        Real  sigma_       ;    /**< @brief 1st gaussian standard deviation @f$ \sigma @f$ (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
        Real  N0_2_        ;    /**< @brief 2nd gaussian @f$ N_0 \left[ m^{-3} \right] @f$. */
        Real  sigma_2_     ;    /**< @brief 2nd gaussian standard deviation @f$ \sigma @f$ (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
        Real  shift_2_     ;    /**< @brief 2nd gaussian shift with respect to the 1st gaussian electric potential (normalized by @f$ -Q @f$) @f$ \left[ V \right] @f$. */
        Real  N0_3_        ;    /**< @brief 3rd gaussian @f$ N_0 \left[ m^{-3} \right] @f$. */
        Real  sigma_3_     ;    /**< @brief 3rd gaussian standard deviation @f$ \sigma @f$ (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
        Real  shift_3_     ;    /**< @brief 3rd gaussian shift with respect to the 1st gaussian electric potential (normalized by @f$ -Q @f$) @f$ \left[ V \right] @f$. */
        Real  N0_4_        ;    /**< @brief 4th gaussian @f$ N_0 \left[ m^{-3} \right] @f$. */
        Real  sigma_4_     ;    /**< @brief 4th gaussian standard deviation @f$ \sigma @f$ (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
        Real  shift_4_     ;    /**< @brief 4th gaussian shift with respect to the 1st gaussian electric potential (normalized by @f$ -Q @f$) @f$ \left[ V \right] @f$. */
        Real  N0_exp_      ;    /**< @brief Exponential @f$ N_0 \left[ m^{-3} \right] @f$. */
        Real  lambda_exp_  ;    /**< @brief Exponential @f$ \lambda @f$ (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
        Real  A_semic_     ;    /**< @brief Area of the semiconductor @f$ \left[ m^2 \right] @f$. */
        Real  C_sb_        ;    /**< @brief Equivalent capacitance of the series semiconductor-insulator (stray-border capacitance) @f$ \left[ F \right] @f$. */
        Index nNodes_      ;    /**< @brief Number of nodes that form the mesh. */
        Index nSteps_      ;    /**< @brief Number of steps to simulate. */
        Real  V_min_       ;    /**< @brief Minimum voltage @f$ \left[ V \right] @f$. */
        Real  V_max_       ;    /**< @brief Maximum voltage @f$ \left[ V \right] @f$. */
};

// Implementations.
inline const Index & ParamList::simulationNo() const
{
    return simulationNo_;
}

inline const Real & ParamList::t_semic() const
{
    return t_semic_;
}

inline const Real & ParamList::t_ins() const
{
    return t_ins_;
}

inline const Real & ParamList::eps_semic() const
{
    return eps_semic_;
}

inline const Real & ParamList::eps_ins() const
{
    return eps_ins_;
}

inline const Real & ParamList::Wf() const
{
    return Wf_;
}

inline const Real & ParamList::Ea() const
{
    return Ea_;
}

inline const Real & ParamList::N0() const
{
    return N0_;
}

inline const Real & ParamList::sigma() const
{
    return sigma_;
}

inline const Real & ParamList::N0_2() const
{
    return N0_2_;
}

inline const Real & ParamList::sigma_2() const
{
    return sigma_2_;
}

inline const Real & ParamList::shift_2() const
{
    return shift_2_;
}

inline const Real & ParamList::N0_3() const
{
    return N0_3_;
}

inline const Real & ParamList::sigma_3() const
{
    return sigma_3_;
}

inline const Real & ParamList::shift_3() const
{
    return shift_3_;
}

inline const Real & ParamList::N0_4() const
{
    return N0_4_;
}

inline const Real & ParamList::sigma_4() const
{
    return sigma_4_;
}

inline const Real & ParamList::shift_4() const
{
    return shift_4_;
}

inline const Real & ParamList::N0_exp() const
{
    return N0_exp_;
}

inline const Real & ParamList::lambda_exp() const
{
    return lambda_exp_;
}

inline const Real & ParamList::A_semic() const
{
    return A_semic_;
}

inline const Real & ParamList::C_sb() const
{
    return C_sb_;
}

inline const Index & ParamList::nNodes() const
{
    return nNodes_;
}

inline const Index & ParamList::nSteps() const
{
    return nSteps_;
}

inline const Real & ParamList::V_min() const
{
    return V_min_;
}

inline const Real & ParamList::V_max() const
{
    return V_max_;
}

#endif /* PARAMLIST_H */
