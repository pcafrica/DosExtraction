#include "typedefs.h"

void utility::print_block(const char * string, std::ostream & os)
{
	os << "╔";
	
	for ( unsigned i = 0; i < strlen(string); ++i ) {
		os << "═";
	}
	
	os << "╗" << std::endl;
	
	os << "║" << string << "║" << std::endl;
	
	os << "╚";
	
	for ( unsigned i = 0; i < strlen(string); ++i ) {
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
