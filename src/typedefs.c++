#include "typedefs.h"

namespace utility
{
	std::string full_path(const std::string & filename, const std::string & relative_directory)
	{
		return (filename[0] == '/') ? filename : (relative_directory + filename);
	}
	
	void print_block(const char * string, std::ostream & os)
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
	
	void print_done(std::ostream & os)
	{
		os << " DONE! ✓" << std::endl;
		
		return;
	}
}
