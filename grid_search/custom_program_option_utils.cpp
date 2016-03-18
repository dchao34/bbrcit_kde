#include <string>
#include <fstream>
#include <algorithm>

#include "custom_program_option_utils.h"

using namespace std;

template<>
std::vector<double> tokenize<double>(const string &line, string delims) {
  vector<string> string_result = tokenize<string>(line, delims);
  vector<double> result;
  for (auto &s : string_result) { result.push_back(std::stod(s)); }
  return result;
}

template<>
vector<string> tokenize<string>(const string &line, string delims) {

  vector<string> result;

  for (auto bit = line.begin();;) {
    auto eit = find_first_of(bit, line.end(), delims.begin(), delims.end());

    // invariant at this line: 
    // 1. bit points to the start of the next token. 
    // 2. eit points to one past the end of the next token. 
    if (bit == eit) {
      result.push_back("");
    } else {
      result.push_back(string(bit, eit));
    }

    if (eit == line.end()) {
      break;
    } else {
      bit = eit + 1;
    }
  }

  return result;
}

