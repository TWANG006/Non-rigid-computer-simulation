
#ifndef _UTILITIES_H
#define	_UTILITIES_H

#include <iostream>		//	standard iostream operations
#include <limits>		//	numeric limits
#include <cmath>		//	mathematical functions
#include <cstdlib>		//	C-string functions
#include <cstddef>		//	C library language support
#include <fstream>		//	file input and ouput
#include <cctype>		//	character classification
#include <ctime>		//	data and time functions
#include <vector>
#include <map>

bool user_say_yes();
enum Error_code
{
	success,fail,ranges_error,underflow,overflow,fatal,
	not_present,cuplicate_error,entry_inserted,entry_found,
	internal_error
};

#endif // !_UTILITIES_H