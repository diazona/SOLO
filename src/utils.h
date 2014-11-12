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
#include <sstream>
#include <string>
#include <vector>

/**
 * Trim characters from a string, or from part of a string.
 * 
 * The effective functionality of `trim is as follows:
 * 
 * - Take the substring of `str` from position `begin` (inclusive) to
 *   `end` (exclusive)
 * - Remove from the beginning of that substring as many of the characters
 *   in `chars` as possible,  until a character not in `chars` is encountered
 * - Remove from the end of the substring in the same way
 * - Replace the substring from the original string with the trimmed substring
 * 
 * except that `str` is not modified; `trim` creates a new `string` object
 * and returns it.
 * 
 * @param str the string to trim
 * @param chars the characters to trim, defaulting to whitespace (`" \r\n\t"`)
 * @param begin the position to start trimming from, defaulting to `0`
 * @param end the position to finish trimming, defaulting to the end of the string
 */
std::string trim(const std::string& str, const std::string& chars = " \r\n\t", size_t begin = 0, size_t end = std::string::npos);

/**
 * Joins strings or string representations of things together, separated by a
 * delimiter.
 */
template<class InputIterator>
std::string join(const std::string& delim, InputIterator begin, InputIterator end) {
    std::ostringstream oss;
    if (begin != end) {
        oss << *begin;
        while (begin != end) {
            oss << delim;
            oss << *begin;
        }
    }
    return oss.str();
}

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
