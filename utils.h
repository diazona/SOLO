#ifndef _UTILS_H_
#define _UTILS_H_

#include <algorithm>
#include <string>
#include <vector>

std::string trim(const std::string& str, const std::string& chars, size_t begin = 0, size_t end = std::string::npos);

std::vector<std::string> split(const std::string& str, const std::string& chars, size_t limit = std::string::npos, size_t i_begin = 0, size_t i_end = std::string::npos);

template<typename T>
T min(std::vector<T> v) {
    return *(std::min_element(v.begin(), v.end()));
}

template<typename T>
T max(std::vector<T> v) {
    return *(std::max_element(v.begin(), v.end()));
}

#endif
