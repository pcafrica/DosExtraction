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
	
	if ( !input_.is_open() ) {
		throw std::ifstream::failure("ERROR: input file cannot be read or wrong filename provided.");
	}
	
	// Get number of rows.
	while ( std::getline(input_, line_, '\n') ) {
		++nRows_;
	}
	
	if ( hasHeaders_ ) {
		--nRows_;
	}
	
	reset();
	
	std::getline(input_, line_, '\n');	// Import first row.
	
	// Check separator.
	if ( std::find(line_.begin(), line_.end(), ',') != line_.end() ) {
		separator_ = ',';
	} else if ( std::find(line_.begin(), line_.end(), '\t') != line_.end() ) {
		separator_ = '\t';
	} else if ( std::find(line_.begin(), line_.end(), ':') != line_.end() ) {
		separator_ = ':';
	} else if ( std::find(line_.begin(), line_.end(), ' ') != line_.end() ) {
		separator_ = ' ';
	} else {
		throw std::ifstream::failure("ERROR: input file isn't either comma-, TAB-, colon- or space-separated.");
	}
	
	// Get number of columns.
	std::stringstream line_stream(line_);	// "line_" contains first row.
	std::string string;
	
	// std::getline(line_stream, string, separator_);	// Skip first column (containing the simulation index).
	
	while ( std::getline(line_stream, string, separator_) ) {
		++nCols_;
	}
	
	reset();
}

RowVectorXd CsvParser::importRow(const unsigned & index)
{
	assert( index >= 1 && index <= nRows_ );
	
	reset();
	
	RowVectorXd data( nCols_ );
	
	for ( unsigned i = 0; i < index; ++i ) {
		std::getline(input_, line_, '\n');	// Read "i"-th row.
	}
	
	// "line_" is now containing "index"-th row.
	
	// Start import.
	{
		std::stringstream line_stream(line_);
		std::string string;
		
		for ( unsigned j = 0; j < nCols_; ++j ) {
			std::getline(line_stream, string, separator_);
			
			data(j) = (double) atof(string.c_str());
		}
	}
	
	return data;
}

MatrixXd CsvParser::importRows(const std::initializer_list<unsigned> & indexes)
{
	assert( indexes.size() > 0 );
	
	MatrixXd Data( indexes.size(), nCols_ );
	
	unsigned i = 0;
	
	for ( unsigned index : indexes ) {
		Data.row(i) = importRow( index );
		++i;
	}
	
	return Data;
}

MatrixXd CsvParser::importFirstRows(const unsigned & nRows)
{
	assert( nRows >= 1 && nRows <= nRows_ );
	
	MatrixXd Data( nRows, nCols_ );
	
	for ( int i = 0; i < Data.rows(); ++i ) {
		Data.row(i) = importRow(i + 1);
	}
	
	return Data;
}

VectorXd CsvParser::importCol(const unsigned & index)
{
	assert( index >= 1 && index <= nCols_ );
	
	reset();
	
	VectorXd data( nRows_ );
	
	// Start import.
	for ( unsigned i = 0; i < nRows_; ++i ) {	// For each row.
		std::getline(input_, line_, '\n');	// Read "i"-th row.
		
		std::stringstream line_stream(line_);	// "line_" already containing first row;
		std::string string;
		
		for ( unsigned j = 0; j < index; ++j ) {
			std::getline(line_stream, string, separator_);
		}
		
		data(i) = (double) atof(string.c_str());
	}
	
	return data;
}

MatrixXd CsvParser::importCols(const std::initializer_list<unsigned> & indexes)
{
	assert( indexes.size() > 0 );
	
	MatrixXd Data( nRows_, indexes.size() );
	
	unsigned j = 0;
	
	for ( int index : indexes ) {
		Data.col(j) = importCol( index );
		++j;
	}
	
	return Data;
}

MatrixXd CsvParser::importFirstCols(const unsigned & nCols)
{
	assert( nCols >= 1 && nCols <= nCols_ );
	
	MatrixXd Data( nRows_, nCols );
	
	for ( int j = 0; j < Data.cols(); ++j ) {
		Data.col(j) = importRow(j + 1);
	}
	
	return Data;
}

double CsvParser::importCell(const unsigned & rowIndex, const unsigned & colIndex)
{
	assert( rowIndex >= 1 && rowIndex <= nRows_ );
	assert( colIndex >= 1 && colIndex <= nCols_ );
	
	return importRow(rowIndex)(colIndex - 1);
}

MatrixXd CsvParser::importAll()
{
	return importFirstRows( nRows_ );
}

void CsvParser::reset()
{
	input_.clear();	// Reset eof flag.
	input_.seekg(0, std::ios::beg);	// Go to beginning of file.
	
	if ( hasHeaders_ ) {	// Ignore row containing headers.
		std::getline(input_, line_, '\n');
	}
	
	return;
}
