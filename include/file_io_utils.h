#ifndef _FILE_IO_UTILS_H_
#define _FILE_IO_UTILS_H_

#include <fstream>
#include <string>
#include <vector>

void open_for_reading(std::ifstream&, std::string);
void open_for_writing(std::ofstream&, std::string);
void tokenize(const std::string &line, 
              std::vector<std::string> &result, 
              std::string delims);

#endif
