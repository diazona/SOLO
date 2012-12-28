#include <cassert>
#include "utils.h"

using namespace std;

// Trim given characters from the beginning and end of a string.
string trim(const string& str, const string& chars, size_t begin, size_t end) {
    size_t l_begin = str.find_first_not_of(chars, begin);
    size_t l_end = str.find_last_not_of(chars, end);
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

int main(int argc, char** argv) {
    cout << "trim(' abcd e ', ' ') = '" << trim(" abcd e ", " ") << "'" << endl;
    cout << "trim('   abcd e   ', ' ') = '" << trim("   abcd e   ", " ") << "'" << endl;
    cout << "trim('abcdefghijk', 'afk') = '" << trim("abcdefghijk", "afk") << "'" << endl;
    
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
    return 0;
}
#endif