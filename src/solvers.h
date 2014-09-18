/* C++ */

/**
 * @file   solvers.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * @brief Generic solvers for PDEs.
 *
 */

#ifndef SOLVERS_H
#define SOLVERS_H

#include "charge.h"
#include "typedefs.h"

#include <utility>	// std::make_pair
#include <limits>	// std::numeric_limits<>::epsilon

class NonLinearPoisson1D;	// Forward declaration.

/**
 * @class PdeSolver1D
 *
 * Matrices are held in a sparse format.
 *
 * @brief Abstract class providing methods to assemble matrices to solve one-dimensional PDEs.
 *
 */
class PdeSolver1D
{
	public:
		friend class NonLinearPoisson1D;
		// Now NonLinearPoisson1D can access system matrices with no need to
		// copy them through getter methods that slow down the program execution.
		
		/**
		 * @brief Default constructor (deleted since it is required to specify the mesh).
		 */
		PdeSolver1D() = delete;
		/**
		 * @brief Constructor.
		 * @param[in] mesh : the mesh.
		 */
		PdeSolver1D(VectorXd &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~PdeSolver1D() = default;
		
		/**
		 * Build the matrix for the advection-diffusion problem:
		 * @f$ -\nabla\cdot\left(\alpha\cdot\gamma\left(\eta\nabla u - \beta u\right)\right) = f @f$.
		 *
		 * @brief Assemble the matrix for an advection-diffusion term.
		 * @param[in] alpha : @f$ \alpha @f$, an element-wise constant function;
		 * @param[in] gamma : @f$ \gamma @f$, an element-wise linear   function;
		 * @param[in] eta   : @f$ \eta   @f$, an element-wise linear   function;
		 * @param[in] beta  : @f$ \beta  @f$, an element-wise constant function.
		 */
		virtual void assembleAdvDiff(const VectorXd & alpha, const VectorXd & gamma,
		                             const VectorXd & eta, const VectorXd & beta) = 0;
		/**
		 * Build the matrix for the diffusion problem:
		 * @f$ -\nabla\cdot\left(\epsilon\cdot\kappa\nabla u\right) = f @f$.
		 *
		 * @brief Assemble the matrix for a diffusion term.
		 * @param[in] eps   : @f$ \epsilon @f$, an element-wise constant function;
		 * @param[in] kappa : @f$ \kappa   @f$, an element-wise linear   function.
		 */
		virtual void assembleStiff  (const VectorXd & eps, const VectorXd & kappa) = 0;
		/**
		 * Build the mass matrix for the reaction problem:
		 * @f$ \delta\cdot\zeta u = f @f$.
		 *
		 * @brief Assemble the matrix for a reaction term.
		 * @param[in] delta : @f$ \delta @f$, an element-wise constant function;
		 * @param[in] zeta  : @f$ \zeta  @f$, an element-wise linear   function.
		 */
		virtual void assembleMass   (const VectorXd & delta, const VectorXd & zeta) = 0;
		
		/**
		 * @name Getter methods
		 * @{
		 */
		inline const SparseXd & AdvDiff() const;
		inline const SparseXd & Stiff  () const;
		inline const SparseXd & Mass   () const;
		/**
		 * @}
		 */
		
	protected:
		VectorXd mesh_  ;	/**< @brief The mesh. */
		unsigned nNodes_;	/**< @brief Number of nodes that form the mesh. */
		
		SparseXd AdvDiff_;	/**< @brief Matrix for an advection-diffusion term. */
		SparseXd Stiff_  ;	/**< @brief Stiffness matrix. */
		SparseXd Mass_   ;	/**< @brief Mass matrix. */
};

/**
 * @class Bim1D
 *
 * Matrices are held in a sparse format.
 *
 * @brief Class derived from @ref PdeSolver1D, providing a finite volume Box Integration Method (BIM) solver.
 *
 */
class Bim1D : public PdeSolver1D
{
	public:
		/**
		 * @brief Default constructor (deleted since it is required to specify the mesh).
		 */
		Bim1D() = delete;
		/**
		 * @brief Constructor.
		 * @param[in] mesh : the mesh coordinates.
		 */
		Bim1D(VectorXd &);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~Bim1D() = default;
		
		/**
		 * @f[ M_{log}(x_1, x_2) = \frac{x_2 - x_1}{\log{x_2} - \log{x_1}} =
		 * \frac{x_2 - x_1}{\log\left(\frac{x_2}{x_1}\right)} ~ .@f]
		 * @brief Compute the element-wise logarithmic mean of two vectors.
		 *
		 * @param[in] x1 : the first vector;
		 * @param[in] x2 : the second vector.
		 * @returns the vector of the logarithmic means.
		 */
		static VectorXd log_mean(const VectorXd &, const VectorXd &);
		/**
		 * @f[ \mathfrak{B}(x) = \frac{x}{e^x - 1} ~ .@f]
		 * @brief Compute the values of the Bernoulli function.
		 * @param[in] x : the vector of the values to compute the Bernoulli function at.
		 * @returns the pair @f$ \left(\mathfrak{B}(x), \mathfrak{B}(-x)\right) @f$.
		 */
		static std::pair<VectorXd, VectorXd> bernoulli(const VectorXd &);
		
		/**
		 * Build the Scharfetter-Gummel stabilized stiffness matrix for:
		 * @f$ -\nabla\cdot\left(\alpha\cdot\gamma\left(\eta\nabla u - \beta u\right)\right) = f @f$.
		 *
		 * @brief Assemble the matrix for an advection-diffusion term.
		 * @param[in] alpha : @f$ \alpha @f$, an element-wise constant function;
		 * @param[in] gamma : @f$ \gamma @f$, an element-wise linear   function;
		 * @param[in] eta   : @f$ \eta   @f$, an element-wise linear   function;
		 * @param[in] beta  : @f$ \beta  @f$, an element-wise constant function.
		 */
		virtual void assembleAdvDiff(const VectorXd &, const VectorXd &, const VectorXd &, const VectorXd &) override;
		/**
		 * Build the standard finite element stiffness matrix for the diffusion problem:
		 * @f$ -\nabla\cdot\left(\epsilon\cdot\kappa\nabla u\right) = f @f$.
		 *
		 * @brief Assemble the matrix for a diffusion term.
		 * @param[in] eps   : @f$ \epsilon @f$, an element-wise constant function;
		 * @param[in] kappa : @f$ \kappa   @f$, an element-wise linear   function.
		 */
		virtual void assembleStiff  (const VectorXd &, const VectorXd &) override;
		/**
		 * Build the lumped finite element mass matrix for the reaction problem:
		 * @f$ \delta\cdot\zeta u = f @f$.
		 *
		 * @brief Assemble the matrix for a reaction term.
		 * @param[in] delta : @f$ \delta @f$, an element-wise constant function;
		 * @param[in] zeta  : @f$ \zeta  @f$, an element-wise linear   function.
		 */
		virtual void assembleMass   (const VectorXd &, const VectorXd &) override;
};

/**
 * @class NonLinearPoisson1D
 *
 * A Newton method is applied in order to solve:
 * @f[ -\frac{\mathrm{d}}{\mathrm{d}z} \left(\epsilon(z) \cdot \frac{\mathrm{d}\varphi}{\mathrm{d}z}(z) \right) =
 * - q \cdot \frac{N_0}{\sqrt{\pi}} \int_{-\infty}^{+\infty} \exp\left(-\alpha^2\right) \left( 1 +
 * \exp\left( \frac{\sqrt{2}\sigma\alpha - q\varphi(z)}{K_B \cdot T} \right) \right)^{-1} \mathrm{d}\alpha ~ .
 * @f]
 *
 *
 * @brief Provide a solver for a non-linear Poisson equation.
 *
 */
class NonLinearPoisson1D
{
	public:
		/**
		 * @brief Default constructor (deleted since it is required to specify the solver to be used).
		 */
		NonLinearPoisson1D() = delete;
		/**
		 * @brief Constructor.
		 * @param[in] solver          : the solver to be used;
		 * @param[in] maxIterationsNo : maximum number of iterations desired;
		 * @param[in] tolerance       : tolerance desired.
		 */
		NonLinearPoisson1D(const PdeSolver1D &, const unsigned & = 100, const double & = 1.0e-6);
		/**
		 * @brief Destructor (defaulted).
		 */
		virtual ~NonLinearPoisson1D() = default;
		
		/**
		 * @brief Apply a Newton method to the equation and then discretize it using the solver specified.
		 * @param[in] mesh       : the mesh;
		 * @param[in] init_guess : initial guess for the Newton algorithm;
		 * @param[in] charge_fun : an object of class @ref Charge specifying how to compute total electric charge.
		 */
		void apply(const VectorXd &, const VectorXd &, Charge &);
		
		/**
		 * @name Getter methods
		 * @{
		 */
		inline const VectorXd & phi()    const;
		inline const VectorXd & norm()   const;
		inline const double   & qTot()   const;
		inline const double   & cTot()   const;
		/* inline const double   & cTot_n() const; */
		/**
		 * @}
		 */
		
	private:
		/**
		 * @brief Compute the Jacobi matrix.
		 * @param[in] x : the vector where to start from.
		 * @returns the Jacobi matrix in a sparse format.
		 */
		SparseXd computeJac(const VectorXd &) const;
		
		const PdeSolver1D & solver_;	/**< @brief Solver handler. */
		
		unsigned maxIterationsNo_;	/**< @brief Maximum number of iterations. */
		double   tolerance_      ;	/**< @brief Tolerance. */
		
		VectorXd phi_ ;	/**< @brief The electric potential. */
		VectorXd norm_;	/**< @brief Vector holding @f$ L^\infty @f$-norm errors for each iteration. */
		
		double qTot_  ;	/**< @brief Total charge. */
		double cTot_  ;	/**< @brief Total capacitance. */
		/* double cTot_n_; */
};

// Implementations.
inline const SparseXd & PdeSolver1D::AdvDiff() const
{
	return AdvDiff_;
}

inline const SparseXd & PdeSolver1D::Stiff() const
{
	return Stiff_;
}

inline const SparseXd & PdeSolver1D::Mass() const
{
	return Mass_;
}

inline const VectorXd & NonLinearPoisson1D::phi() const
{
	return phi_;
}

inline const VectorXd & NonLinearPoisson1D::norm() const
{
	return norm_;
}

inline const double & NonLinearPoisson1D::qTot() const
{
	return qTot_;
}

inline const double & NonLinearPoisson1D::cTot() const
{
	return cTot_;
}

/* inline const double & NonLinearPoisson1D::cTot_n() const
{
	return cTot_n_;
} */

#endif /* SOLVERS_H */
