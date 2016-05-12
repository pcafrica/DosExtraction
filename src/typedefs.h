/* C++11 */

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

#include "GetPot.h"

#include <iostream>
#include <fstream>

typedef double Real;    /**< @brief Typedef for real numbers. */
typedef ptrdiff_t Index;    /**< @brief Typedef for indexing variables. */

#include "physicalConstants.h"

using namespace Eigen;

/**
 * @brief Template alias for @ref Eigen matrices.
 * @tparam ScalarType : the scalar type.
 */
template<typename ScalarType>
using MatrixX = Matrix<ScalarType, Dynamic, Dynamic>;

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

using MatrixXr    = MatrixX<Real>            ;    /**< @brief Typedef for dense real-valued dynamic-sized matrices. */
using VectorXr    = VectorX<Real>            ;    /**< @brief Typedef for dense real-valued dynamic-sized column vectors. */
using RowVectorXr = Matrix<Real, 1, Dynamic> ;    /**< @brief Typedef for dense real-valued dynamic-sized row vectors. */

using SparseXr = SparseMatrix<Real>;    /**< @brief Typedef for sparse real-valued dynamic-sized matrices. */

namespace constants
{
    const Index PARAMS_NO = 27;    /**< @brief Number of parameters required in input file. */
    
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
    
    /**
     * @brief Auxiliary function to write an Eigen matrix to a binary file.
     * @tparam Matrix : the matrix type.
     * @param[in] filename : the filename;
     * @param[in] matrix   : the matrix to write.
     */
    template<class Matrix>
    void write_binary(const std::string &, const Matrix &);
    
    /**
     * @brief Auxiliary function to write a binary file to an Eigen matrix.
     * @tparam Matrix : the matrix type.
     * @param[in] filename : the filename;
     * @param[in] matrix   : the matrix to write on.
     */
    template<class Matrix>
    void read_binary(const std::string &, Matrix &);
}

// Implementations.
template<class Matrix>
void utility::write_binary(const std::string & filename, const Matrix & matrix)
{
    std::ofstream output;
    
    output.open(filename, std::ios_base::out | std::ios_base::binary);
    
    if (output.bad())
    {
        throw std::ofstream::failure ("ERROR: output files cannot be opened or directory does not exist.");
    }
    
    typename Matrix::Index nRows = matrix.rows();
    typename Matrix::Index nCols = matrix.cols();
    
    output.write((char*) (&nRows), sizeof(typename Matrix::Index));
    output.write((char*) (&nCols), sizeof(typename Matrix::Index));
    output.write((char*) matrix.data(), nRows * nCols * sizeof(typename Matrix::Scalar));
    
    output.close();
}

template<class Matrix>
void utility::read_binary(const std::string & filename, Matrix & matrix)
{
    std::ifstream input;
    
    input.open(filename, std::ios_base::in | std::ios_base::binary);
    
    if (input.bad())
    {
        throw std::ofstream::failure ("ERROR: input files cannot be opened or directory does not exist.");
    }
    
    typename Matrix::Index nRows = 0, nCols = 0;
    
    input.read((char*) (&nRows), sizeof(typename Matrix::Index));
    input.read((char*) (&nCols), sizeof(typename Matrix::Index));
    
    matrix.resize(nRows, nCols);
    
    input.read((char*) matrix.data() , nRows * nCols * sizeof(typename Matrix::Scalar));
    
    input.close();
}

#endif /* TYPEDEFS_H */
