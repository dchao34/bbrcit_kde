#ifndef _CUSTOM_PROGRAM_OPTION_UTILS_H_
#define _CUSTOM_PROGRAM_OPTION_UTILS_H_

#include <ostream>
#include <vector>
#include <algorithm>

// given a string (`line`) delimited with some set of characters (`delims`), 
// return a vector storing each token in the order they are encountered.
template<typename T> 
std::vector<T> tokenize(const std::string &line, std::string delims=" ");

template<> 
std::vector<std::string> 
tokenize<std::string>(const std::string &line, std::string delims);

template<> 
std::vector<double> 
tokenize<double>(const std::string &line, std::string delims);

#endif
