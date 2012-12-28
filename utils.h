#ifndef _UTILS_H_
#define _UTILS_H_

#include <string>
#include <vector>

std::string trim(const std::string& str, const std::string& chars, size_t begin = 0, size_t end = std::string::npos);

std::vector<std::string> split(const std::string& str, const std::string& chars, size_t limit = std::string::npos, size_t i_begin = 0, size_t i_end = std::string::npos);

#endif
