/* C++ */

/**
 * @file   paramList.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief List of simulation parameters.
 *
 */

#ifndef PARAMLIST_H
#define PARAMLIST_H

#include "typedefs.h"

/**
 * @class ParamList
 *
 * It can include up to 4 gaussians, later combined to compute total charge.
 *
 * @brief Class providing methods to handle a list of parameters.
 *
 */
class ParamList
{
	public:
		friend class GaussianCharge;
		friend class       DosModel;
		// Now GaussianCharge and DosModel can access private parameters with no need to
		// copy them through getter methods that slow down the program execution.
		
		/**
		 * @brief Default constructor (defaulted).
		 */
		ParamList() = default;
		/**
		 * @brief Explicit conversion constructor.
		 * @param[in] list : a row vector containing a parameters list (for example got by a @ref CsvParser object).
		 * Parameters should be sorted in the same order as specified above.
		 */
		explicit ParamList(const RowVectorXd &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~ParamList() = default;
		
		/**
		 * @name Getter methods.
		 * @{
		 */
		inline const unsigned & simulationNo() const;
		inline const double   & t_semic()      const;
		inline const double   & t_ins()        const;
		inline const double   & eps_semic()    const;
		inline const double   & eps_ins()      const;
		inline const double   & Wf()           const;
		inline const double   & Ea()           const;
		inline const double   & N0()           const;
		inline const double   & sigma()        const;
		inline const double   & N0_2()         const;
		inline const double   & sigma_2()      const;
		inline const double   & shift_2()      const;
		inline const double   & N0_3()         const;
		inline const double   & sigma_3()      const;
		inline const double   & shift_3()      const;
		inline const double   & N0_4()         const;
		inline const double   & sigma_4()      const;
		inline const double   & shift_4()      const;
		inline const unsigned & nNodes()       const;
		inline const unsigned & nSteps()       const;
		inline const double   & V_min()        const;
		inline const double   & V_max()        const;
		/**
		 * @}
		 */
		
	private:
		unsigned simulationNo_;	/**< @brief Index of the simulation. */
		double   t_semic_     ;	/**< @brief Semiconductor layer thickness @f$ \left[ m \right] @f$. */
		double   t_ins_       ;	/**< @brief Insulator layer thickness @f$ \left[ m \right] @f$. */
		double   eps_semic_   ;	/**< @brief Semiconductor layer relative electrical permittivity @f$ \left[ ~ \right] @f$. */
		double   eps_ins_     ;	/**< @brief Insulator layer relative electrical permittivity @f$ \left[ ~ \right] @f$. */
		double   Wf_          ;	/**< @brief Work-function @f$ \left[ V \right] @f$. */
		double   Ea_          ;	/**< @brief Electron affinity @f$ \left[ V \right] @f$. */
		double   N0_          ;	/**< @brief 1st gaussian mean @f$ \left[ m^{-3} \right] @f$. */
		double   sigma_       ;	/**< @brief 1st gaussian standard deviation (normalized by @f$ K_B \cdot T @f$) @f$ \left[ ~ \right] @f$. */
		double   N0_2_        ;	/**< @brief 2nd gaussian mean. */
		double   sigma_2_     ;	/**< @brief 2nd gaussian standard deviation. */
		double   shift_2_     ;	/**< @brief 2nd gaussian shift with respect to the 1st gaussian electric potential. */
		double   N0_3_        ;	/**< @brief 3rd gaussian mean. */
		double   sigma_3_     ;	/**< @brief 3rd gaussian standard deviation. */
		double   shift_3_     ;	/**< @brief 3rd gaussian shift with respect to the 1st gaussian electric potential. */
		double   N0_4_        ;	/**< @brief 4th gaussian mean. */
		double   sigma_4_     ;	/**< @brief 4th gaussian standard deviation. */
		double   shift_4_     ;	/**< @brief 4th gaussian shift with respect to the 1st gaussian electric potential. */
		unsigned nNodes_      ;	/**< @brief No. of nodes that form the mesh. */
		unsigned nSteps_      ;	/**< @brief No. of steps to simulate. */
		double   V_min_       ;	/**< @brief Minimum voltage @f$ \left[ V \right] @f$. */
		double   V_max_       ;	/**< @brief Maximum voltage @f$ \left[ V \right] @f$. */
};

inline const unsigned & ParamList::simulationNo() const
{
	return simulationNo_;
}

inline const double & ParamList::t_semic() const
{
	return t_semic_;
}

inline const double & ParamList::t_ins() const
{
	return t_ins_;
}

inline const double & ParamList::eps_semic() const
{
	return eps_semic_;
}

inline const double & ParamList::eps_ins() const
{
	return eps_ins_;
}

inline const double & ParamList::Wf() const
{
	return Wf_;
}

inline const double & ParamList::Ea() const
{
	return Ea_;
}

inline const double & ParamList::N0() const
{
	return N0_;
}

inline const double & ParamList::sigma() const
{
	return sigma_;
}

inline const double & ParamList::N0_2() const
{
	return N0_2_;
}

inline const double & ParamList::sigma_2() const
{
	return sigma_2_;
}

inline const double & ParamList::shift_2() const
{
	return shift_2_;
}

inline const double & ParamList::N0_3() const
{
	return N0_3_;
}

inline const double & ParamList::sigma_3() const
{
	return sigma_3_;
}

inline const double & ParamList::shift_3() const
{
	return shift_3_;
}

inline const double & ParamList::N0_4() const
{
	return N0_4_;
}

inline const double & ParamList::sigma_4() const
{
	return sigma_4_;
}

inline const double & ParamList::shift_4() const
{
	return shift_4_;
}

inline const unsigned & ParamList::nNodes() const
{
	return nNodes_;
}

inline const unsigned & ParamList::nSteps() const
{
	return nSteps_;
}

inline const double & ParamList::V_min() const
{
	return V_min_;
}

inline const double & ParamList::V_max() const
{
	return V_max_;
}

#endif /* PARAMLIST_H */
