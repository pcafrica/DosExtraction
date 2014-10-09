/* C++ */

/**
 * @file   csvParser.h
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Tools to store content from a .csv file in matrices or vectors.
 *
 */

#ifndef CSVPARSER_H
#define CSVPARSER_H

#include "typedefs.h"

#include <fstream>
#include <sstream>
#include <string>
#include <utility>    // std::initializer_list<>

/**
 * @class CsvParser
 *
 * @brief Class providing methods to read @b numeric content from a .csv file
 * and to store it in @ref Eigen matrices or vectors.
 *
 */
class CsvParser
{
  public:
    /**
     * @brief Default constructor (deleted since it is required to specify at least a filename).
     */
    CsvParser() = delete;
    /**
     * @brief Constructor: load the input file and check its compatibility with the code.
     * @param[in] input_filename : the name of the input file;
     * @param[in] hasHeaders     : bool to specify if first row contains headers or not;
     * if @b true, first row is always ignored.
     */
    CsvParser(const std::string &, const bool & = true);
    /**
     * @brief Destructor: close the input file.
     */
    virtual ~CsvParser();
    
    /**
     * @name Getter methods
     * @{
     */
    inline const Index & nRows() const;
    inline const Index & nCols() const;
    
    /**
     * @}
     */
    
    /**
     * @brief Method to import a row from the input file.
     * @param[in] index : the row index.
     * @returns a row vector containing the content read.
     */
    RowVectorXr importRow      (const Index                        &);
    /**
     * @brief Method to import multiple rows from the input file.
     * @param[in] indexes : initializer list containing the row indexes (e.g. something like {1, 3, 4}).
     * @returns a matrix containing the content read (row by row).
     */
    MatrixXr    importRows     (const std::initializer_list<Index> &);
    /**
     * @brief Method to import the first @a nRows rows from the input file.
     * @param[in] nRows : the number of rows to import.
     * @returns a matrix containing the content read (row by row).
     */
    MatrixXr    importFirstRows(const Index                        &);
    
    /**
     * @brief Method to import a column from the input file.
     * @param[in] index : the column index.
     * @returns a column vector containing the content read.
     */
    VectorXr    importCol      (const Index                        &);
    /**
     * @brief Method to import multiple columns from the input file.
     * @param[in] indexes : initializer list containing the column indexes (e.g. something like {1, 3, 4}).
     * @returns a matrix containing the content read (column by column).
     */
    MatrixXr    importCols     (const std::initializer_list<Index> &);
    /**
     * @brief Method to import the first @a nCols columns from the input file.
     * @param[in] nCols : the number of columns to import.
     * @returns a matrix containing the content read (column by column).
     */
    MatrixXr    importFirstCols(const Index                        &);
    
    /**
     * @brief Method to import a single cell from the input file.
     * @param[in] rowIndex : the cell row index.
     * @param[in] colIndex : the cell column index.
     * @returns a scalar containing the value read.
     */
    Real      importCell     (const Index &, const Index &);
    
    /**
     * @brief Method to import the whole input file.
     * @returns a matrix containing the content read (cell by cell).
     */
    MatrixXr    importAll      ();
    
  private:
    /**
     * @brief Reset all the flags for @a input_ and go back to the beginning of file (possibly by ignoring headers).
     */
    void reset();
    
    bool  hasHeaders_;    /**< @brief bool to determine if first row contains headers or not. */
    Index nRows_     ;    /**< @brief Number of rows in the input file. */
    Index nCols_     ;    /**< @brief Number of columns in the input file. */
    
    std::ifstream input_;    /**< @brief Input stream to @a input_filename. */
    std::string   line_ ;    /**< @brief Auxiliary variable to store currently processed line. */
    
    char separator_;    /**< @brief The separator character detected. */
};

// Implementations.
inline const Index & CsvParser::nRows() const
{
  return nRows_;
}

inline const Index & CsvParser::nCols() const
{
  return nCols_;
}

#endif /* CSVPARSER_H */
