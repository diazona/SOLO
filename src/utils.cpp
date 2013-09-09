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

#include <cassert>
#include "utils.h"

using namespace std;

// Trim given characters from the beginning and end of a string.
string trim(const string& str, const string& chars, size_t begin, size_t end) {
    size_t l_begin = str.find_first_not_of(chars, begin);
    size_t l_end = str.find_last_not_of(chars, end - 1);
    return str.substr(l_begin, l_end - l_begin + 1);
}

// Split a string using any of the given chars as delimiters. Empty fields are dropped.
vector<string> split(const string& str, const string& chars, size_t limit, size_t i_begin, size_t i_end) {
    vector<string> tokens;
    size_t delimiter, n = 0;
    if (str.size() < i_end) {
        i_end = str.size();
    }
    while (n++ < limit && i_begin < i_end) {
        delimiter = str.find_first_of(chars, i_begin);
        if (delimiter != i_begin) {
            assert(delimiter > i_begin);
            tokens.push_back(str.substr(i_begin, delimiter - i_begin));
        }
        i_begin = str.find_first_not_of(chars, delimiter);
        assert(i_begin == string::npos || i_begin > delimiter);
    }
    if (i_begin < i_end) {
        tokens.push_back(str.substr(i_begin, i_end - i_begin));
    }
    return tokens;
}

#ifdef UTILS_TEST
#include <iostream>

/**
 * Tests some of the functions in this file.
 */
int main(int argc, char** argv) {
    cout << "trim(' abcd e ', ' ') = '" << trim(" abcd e ", " ") << "'" << endl;
    cout << "trim('   abcd e   ', ' ') = '" << trim("   abcd e   ", " ") << "'" << endl;
    cout << "trim('abcdefghijk', 'afk') = '" << trim("abcdefghijk", "afk") << "'" << endl;

    cout << "trim('abcdefghijk', 'f', 0, 6) = '" << trim("abcdefghijk", "f", 0, 6) << "'" << endl;
    cout << "trim('abcdefghijk', 'afk', 0, 6) = '" << trim("abcdefghijk", "afk", 0, 6) << "'" << endl;
    cout << "trim('abcdefghijk', 'ceg', 0, 6) = '" << trim("abcdefghijk", "ceg", 0, 6) << "'" << endl;
    
    vector<string> tokens = split("abc,def,ghi,jkl", ",");
    cout << "split('abc,def,ghi,jkl', ',') = ";
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
        cout << "'" << *it << "' ";
    }
    cout << endl;

    tokens = split("abcdefghijklmnop", "abcfnp");
    cout << "split('abcdefghijklmnop', 'abcfnp') = ";
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
        cout << "'" << *it << "' ";
    }
    cout << endl;

    tokens = split("abcdefghijklmnop", "acefnp", 2);
    cout << "split('abcdefghijklmnop', 'acefnp', 2) = ";
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
        cout << "'" << *it << "' ";
    }
    cout << endl;

    tokens = split("abcdefghijklmnop", "acefnp", 3);
    cout << "split('abcdefghijklmnop', 'acefnp', 3) = ";
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
        cout << "'" << *it << "' ";
    }
    cout << endl;

    tokens = split("abcdefghijklmnop", "acefnp", 7);
    cout << "split('abcdefghijklmnop', 'acefnp', 7) = ";
    for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++) {
        cout << "'" << *it << "' ";
    }
    cout << endl;
    
    vector<double> v;
    for (double d = 0; d < 10; d++) { v.push_back(d); }
    random_shuffle(v.begin(), v.end());
    cout << "max(0,...,10) = " << max(v) << endl;
    cout << "min(0,...,10) = " << min(v) << endl;
    
    return 0;
}
#endif