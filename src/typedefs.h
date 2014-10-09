/* C++ */

/**
 * @file   typedefs.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Typedefs and utility functions.
 *
 */

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "GetPot"

#include <iostream>
#include <fstream>

#define Real  double    /**< @brief Pre-processor macro for real numbers. */
#define Index ptrdiff_t    /**< @brief Pre-processor macro for indexing variables. */

#include "physicalConstants.h"

using namespace Eigen;

typedef Matrix<Real, Dynamic, Dynamic> MatrixXr   ;    /**< @brief Typedef for dense real-valued dynamic-sized matrices. */
typedef Matrix<Real, Dynamic,       1> VectorXr   ;    /**< @brief Typedef for dense real-valued dynamic-sized column vectors. */
typedef Matrix<Real,       1, Dynamic> RowVectorXr;    /**< @brief Typedef for dense real-valued dynamic-sized row vectors. */

typedef SparseMatrix<Real>             SparseXr   ;    /**< @brief Typedef for sparse real-valued dynamic-sized matrices. */

/**
 * @brief Template alias for @ref Eigen vectors.
 * @tparam ScalarType : the scalar type.
 */
template<typename ScalarType>
using VectorX = Matrix<ScalarType, Dynamic, 1>;

/**
 * @brief Template alias for an @ref Eigen vector of pairs: (@a ScalarType, @a Index).
 * @tparam ScalarType : the scalar type.
 */
template<typename T>
using VectorXpair = VectorX<std::pair<T, Index> >;

namespace constants
{
  const Index PARAMS_NO = 26;    /**< @brief Number of parameters required in input file. */
  
  // From <cmath> library.
  const Real      PI = M_PI              ;    /**< @brief @f$ \pi @f$. */
  const Real SQRT_PI = std::sqrt(PI)     ;    /**< @brief @f$ \sqrt{\pi} @f$. */
  const Real   PI_M4 = 0.7511255444649425;    /**< @brief @f$ \pi^{-\frac{1}{4}} @f$. */
  const Real  SQRT_2 = std::sqrt(2)      ;    /**< @brief @f$ \sqrt{2} @f$. */
}

/**
 * @namespace utility
 *
 * @brief Namespace for utilities and auxiliary functions.
 *
 */
namespace utility
{
  /**
   * @brief Auxiliary function to return the full path to a file.
   * @param[in] filename           : the filename;
   * @param[in] relative_directory : the directory for a relative path.
   * @returns the variable @a filename, if it contains an absolute path; otherwise returns the concatenation
   * of @a relative_directory and @a filename (i.e. the relative path to @a filename).
   */
  std::string full_path(const std::string &, const std::string &);
  
  /**
   * @brief Auxiliary function to print a string inside a block.
   * @param[in]  string : the string to print;
   * @param[out] os     : output stream.
   */
  void print_block(const char *, std::ostream & = std::cout);
  /**
   * @brief Auxiliary function to print a "DONE!" string.
   * @param[out] os : output stream.
   */
  void print_done (std::ostream & = std::cout);
}

#endif /* TYPEDEFS_H */
