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

#include <utility>    // std::make_pair
#include <limits>    // std::numeric_limits<>::epsilon

class NonLinearPoisson1D;    // Forward declaration.

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
    PdeSolver1D(VectorXr &);
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
    virtual void assembleAdvDiff(const VectorXr & alpha, const VectorXr & gamma,
                                 const VectorXr & eta, const VectorXr & beta)  = 0;
    /**
     * Build the matrix for the diffusion problem:
     * @f$ -\nabla\cdot\left(\epsilon\cdot\kappa\nabla u\right) = f @f$.
     *
     * @brief Assemble the matrix for a diffusion term.
     * @param[in] eps   : @f$ \epsilon @f$, an element-wise constant function;
     * @param[in] kappa : @f$ \kappa   @f$, an element-wise linear   function.
     */
    virtual void assembleStiff  (const VectorXr & eps, const VectorXr & kappa) = 0;
    /**
     * Build the mass matrix for the reaction problem:
     * @f$ \delta\cdot\zeta u = f @f$.
     *
     * @brief Assemble the matrix for a reaction term.
     * @param[in] delta : @f$ \delta @f$, an element-wise constant function;
     * @param[in] zeta  : @f$ \zeta  @f$, an element-wise linear   function.
     */
    virtual void assembleMass   (const VectorXr & delta, const VectorXr & zeta) = 0;
    
    /**
     * @name Getter methods
     * @{
     */
    inline const SparseXr & AdvDiff() const;
    inline const SparseXr & Stiff  () const;
    inline const SparseXr & Mass   () const;
    
    /**
     * @}
     */
    
  protected:
    VectorXr mesh_  ;    /**< @brief The mesh. */
    Index    nNodes_;    /**< @brief Number of nodes that form the mesh. */
    
    SparseXr AdvDiff_;    /**< @brief Matrix for an advection-diffusion term. */
    SparseXr Stiff_  ;    /**< @brief Stiffness matrix. */
    SparseXr Mass_   ;    /**< @brief Mass matrix. */
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
    Bim1D(VectorXr &);
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
    static VectorXr log_mean(const VectorXr &, const VectorXr &);
    /**
     * @f[ \mathfrak{B}(x) = \frac{x}{e^x - 1} ~ .@f]
     * @brief Compute the values of the Bernoulli function.
     * @param[in] x : the vector of the values to compute the Bernoulli function at.
     * @returns the pair @f$ \left(\mathfrak{B}(x), \mathfrak{B}(-x)\right) @f$.
     */
    static std::pair<VectorXr, VectorXr> bernoulli(const VectorXr &);
    
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
    virtual void assembleAdvDiff(const VectorXr &, const VectorXr &, const VectorXr &, const VectorXr &) override;
    /**
     * Build the standard finite element stiffness matrix for the diffusion problem:
     * @f$ -\nabla\cdot\left(\epsilon\cdot\kappa\nabla u\right) = f @f$.
     *
     * @brief Assemble the matrix for a diffusion term.
     * @param[in] eps   : @f$ \epsilon @f$, an element-wise constant function;
     * @param[in] kappa : @f$ \kappa   @f$, an element-wise linear   function.
     */
    virtual void assembleStiff  (const VectorXr &, const VectorXr &) override;
    /**
     * Build the lumped finite element mass matrix for the reaction problem:
     * @f$ \delta\cdot\zeta u = f @f$.
     *
     * @brief Assemble the matrix for a reaction term.
     * @param[in] delta : @f$ \delta @f$, an element-wise constant function;
     * @param[in] zeta  : @f$ \zeta  @f$, an element-wise linear   function.
     */
    virtual void assembleMass   (const VectorXr &, const VectorXr &) override;
};

/**
 * @class NonLinearPoisson1D
 *
 * A Newton method is applied in order to solve:
 * @f[ -\frac{\mathrm{d}}{\mathrm{d}z} \left(\epsilon(z) \cdot \frac{\mathrm{d}\varphi}{\mathrm{d}z}(z) \right) =
 * - q \cdot \frac{N_0}{\sqrt{\pi}} \int_{-\infty}^{+\infty} \exp\left(-\alpha^2\right) \left( 1 +
 * \exp\left( \frac{\sqrt{2}\sigma\alpha - q\varphi(z)}{K_B \cdot T} \right) \right)^{-1} \mathrm{d}\alpha ~ . @f]
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
    NonLinearPoisson1D(const PdeSolver1D &, const Index & = 100, const Real & = 1.0e-6);
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
    void apply(const VectorXr &, const VectorXr &, Charge &);
    
    /**
     * @name Getter methods
     * @{
     */
    inline const VectorXr & phi()    const;
    inline const VectorXr & norm()   const;
    inline const Real     & qTot()   const;
    inline const Real     & cTot()   const;
    /* inline const Real     & cTot_n() const; */
    
    /**
     * @}
     */
    
  private:
    /**
     * @brief Compute the Jacobi matrix.
     * @param[in] x : the vector where to start from.
     * @returns the Jacobi matrix in a sparse format.
     */
    SparseXr computeJac(const VectorXr &) const;
    
    const PdeSolver1D & solver_;    /**< @brief Solver handler. */
    
    Index maxIterationsNo_;    /**< @brief Maximum number of iterations. */
    Real  tolerance_      ;    /**< @brief Tolerance. */
    
    VectorXr phi_ ;    /**< @brief The electric potential. */
    VectorXr norm_;    /**< @brief Vector holding @f$ L^\infty @f$-norm errors for each iteration. */
    
    Real qTot_  ;    /**< @brief Total charge. */
    Real cTot_  ;    /**< @brief Total capacitance. */
    /* Real cTot_n_; */
};

// Implementations.
inline const SparseXr & PdeSolver1D::AdvDiff() const
{
  return AdvDiff_;
}

inline const SparseXr & PdeSolver1D::Stiff() const
{
  return Stiff_;
}

inline const SparseXr & PdeSolver1D::Mass() const
{
  return Mass_;
}

inline const VectorXr & NonLinearPoisson1D::phi() const
{
  return phi_;
}

inline const VectorXr & NonLinearPoisson1D::norm() const
{
  return norm_;
}

inline const Real & NonLinearPoisson1D::qTot() const
{
  return qTot_;
}

inline const Real & NonLinearPoisson1D::cTot() const
{
  return cTot_;
}

/* inline const Real & NonLinearPoisson1D::cTot_n() const
{
    return cTot_n_;
} */

#endif /* SOLVERS_H */
