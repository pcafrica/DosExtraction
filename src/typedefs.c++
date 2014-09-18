#include "typedefs.h"

std::string utility::full_path(const std::string & filename, const std::string & relative_directory)
{
	return (filename[0] == '/') ? filename : (relative_directory + filename);
}

void utility::print_block(const char * string, std::ostream & os)
{
	os << "╔";
	
	for ( std::size_t i = 0; i < strlen(string); ++i ) {
		os << "═";
	}
	
	os << "╗" << std::endl;
	
	os << "║" << string << "║" << std::endl;
	
	os << "╚";
	
	for ( std::size_t i = 0; i < strlen(string); ++i ) {
		os << "═";
	}
	
	os << "╝" << std::endl;
	
	return;
}

void utility::print_done(std::ostream & os)
{
	os << " DONE! ✓" << std::endl;
	
	return;
}
