/* C++11 */

/**
 * @file   csvParser.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "csvParser.h"

using namespace constants;

CsvParser::~CsvParser()
{
  input_.close();
}

CsvParser::CsvParser(const std::string & input_filename, const bool & hasHeaders)
  : hasHeaders_(hasHeaders), nRows_(0), nCols_(0)
{
  input_.open(input_filename, std::ios::in);
  
  if ( !input_.is_open() )
    {
      throw std::ifstream::failure("ERROR: input file cannot be read or wrong filename provided.");
    }
    
  // Get number of rows.
  
  if ( hasHeaders_ )    // Skip the row containing headers.
    {
      std::getline(input_, line_, '\n');
    }
    
  while ( std::getline(input_, line_, '\n') )
    {
      ++nRows_;
    }
    
  reset();
  
  std::getline(input_, line_, '\n');    // Import first row.
  
  // Check separator.
  if ( std::find(line_.begin(), line_.end(), ',') != line_.end() )
    {
      separator_ = ',';
    }
  else if ( std::find(line_.begin(), line_.end(), '\t') != line_.end() )
    {
      separator_ = '\t';
    }
  else if ( std::find(line_.begin(), line_.end(), ':') != line_.end() )
    {
      separator_ = ':';
    }
  else if ( std::find(line_.begin(), line_.end(), ' ') != line_.end() )
    {
      separator_ = ' ';
    }
  else
    {
      throw std::ifstream::failure("ERROR: input file isn't either comma-, TAB-, colon- or space-separated.");
    }
    
  // Get number of columns.
  {
    // Import the first row: .csv files are formatted so that
    // each row contains the same number of columns.
    std::istringstream line_stream(line_);
    
    std::string field;
    
    while ( std::getline(line_stream, field, separator_) )
      {
        ++nCols_;
      }
      
    // If field is empty, the last column has not been counted,
    // because the end of line has been reached.
    if ( field.empty() )
      {
        ++nCols_;
      }
  }
  
  reset();
}

RowVectorXr CsvParser::importRow(const Index & index)
{
  assert( index >= 1 && index <= nRows_ );
  
  reset();
  
  std::locale locale;
  
  RowVectorXr data = RowVectorXr::Zero( nCols_ );
  
  for ( Index i = 0; i < index; ++i )
    {
      std::getline(input_, line_, '\n');    // Read "i"-th row.
    }
    
  // "line_" is now containing "index"-th row.
  
  // Start import.
  {
    std::istringstream line_stream(line_);
    std::string field;
    
    for ( Index j = 0; j < nCols_; ++j )
      {
        std::getline(line_stream, field, separator_);
        
        data(j) = (Real) std::atof(field.c_str());
      }
  }
  
  return data;
}

MatrixXr CsvParser::importRows(const std::initializer_list<Index> & indexes)
{
  assert( indexes.size() > 0 );
  
  MatrixXr Data = MatrixXr::Zero( indexes.size(), nCols_ );
  
  Index i = 0;
  
  for ( Index index : indexes )
    {
      Data.row(i) = importRow( index );
      ++i;
    }
    
  return Data;
}

MatrixXr CsvParser::importFirstRows(const Index & nRows)
{
  assert( nRows >= 1 && nRows <= nRows_ );
  
  MatrixXr Data = MatrixXr::Zero( nRows, nCols_ );
  
  for ( Index i = 0; i < Data.rows(); ++i )
    {
      Data.row(i) = importRow(i + 1);
    }
    
  return Data;
}

VectorXr CsvParser::importCol(const Index & index)
{
  assert( index >= 1 && index <= nCols_ );
  
  reset();
  
  std::locale locale;
  
  VectorXr data = VectorXr::Zero( nRows_ );
  
  // Start import.
  for ( Index i = 0; i < nRows_; ++i )    // For each row.
    {
      std::getline(input_, line_, '\n');    // Read "i"-th row.
      
      std::istringstream line_stream(line_);
      std::string field;
      
      for ( Index j = 0; j < index; ++j )
        {
          std::getline(line_stream, field, separator_);
        }
        
      data(i) = (Real) std::atof(field.c_str());
    }
    
  return data;
}

MatrixXr CsvParser::importCols(const std::initializer_list<Index> & indexes)
{
  assert( indexes.size() > 0 );
  
  MatrixXr Data = MatrixXr::Zero( nRows_, indexes.size() );
  
  Index j = 0;
  
  for ( Index index : indexes )
    {
      Data.col(j) = importCol( index );
      ++j;
    }
    
  return Data;
}

MatrixXr CsvParser::importFirstCols(const Index & nCols)
{
  assert( nCols >= 1 && nCols <= nCols_ );
  
  MatrixXr Data = MatrixXr::Zero( nRows_, nCols );
  
  for ( Index j = 0; j < Data.cols(); ++j )
    {
      Data.col(j) = importRow(j + 1);
    }
    
  return Data;
}

Real CsvParser::importCell(const Index & rowIndex, const Index & colIndex)
{
  assert( rowIndex >= 1 && rowIndex <= nRows_ );
  assert( colIndex >= 1 && colIndex <= nCols_ );
  
  return importRow(rowIndex)(colIndex - 1);
}

MatrixXr CsvParser::importAll()
{
  return importFirstRows( nRows_ );
}

void CsvParser::reset()
{
  input_.clear();    // Reset eof flag.
  input_.seekg(0, std::ios::beg);    // Go to beginning of file.
  
  if ( hasHeaders_ )    // Skip the row containing headers.
    {
      std::getline(input_, line_, '\n');
    }
    
  return;
}
