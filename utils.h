/*
 * Part of oneloopcalc
 * 
 * Copyright 2012 David Zaslavsky
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _UTILS_H_
#define _UTILS_H_

#include <algorithm>
#include <string>
#include <vector>

/**
 * Trim any of the characters given in `chars` from the beginning and end of the
 * string given in `str` or its substring from `begin` to `end`), and returns
 * the modified string.
 */
std::string trim(const std::string& str, const std::string& chars, size_t begin = 0, size_t end = std::string::npos);

/**
 * Split the string given in `str` (or its substring from `begin` to `end`) using any
 * of the characters given in `chars` as delimiters, but producing no more than
 * `limit` strings total.
 */
std::vector<std::string> split(const std::string& str, const std::string& chars, size_t limit = std::string::npos, size_t i_begin = 0, size_t i_end = std::string::npos);

/**
 * Returns the smallest element of a list.
 */
template<typename T>
T min(std::vector<T> v) {
    return *(std::min_element(v.begin(), v.end()));
}

/**
 * Returns the largest element of a list.
 */
template<typename T>
T max(std::vector<T> v) {
    return *(std::max_element(v.begin(), v.end()));
}

#endif
