#include <cassert>
#include <fstream>
#include <iostream>
#include "configuration.h"
#include "../utils/utils.h"

using std::ios_base;
using std::endl;
using std::istream;
using std::ifstream;
using std::ostream;
using std::multimap;
using std::pair;
using std::string;
using std::vector;

string identity(const string& key) {
    return string(key);
}

Configuration::Configuration(string (*canonicalize)(const string&)) : canonicalize(canonicalize) {
    if (this->canonicalize == NULL) {
        this->canonicalize = identity;
    }
}

Configuration::Configuration(const string& filename, string (*canonicalize)(const string&)) : canonicalize(canonicalize) {
    ifstream in(filename.c_str());
    if (!in.good()) {
        throw ios_base::failure("Unable to read file " + filename);
    }
    read_config(in);
    in.close();
    if (this->canonicalize == NULL) {
        this->canonicalize = identity;
    }
}

typedef multimap<string,string>::const_iterator option_iterator_type;

string Configuration::get(const string& key, const size_t index) const {
    pair<option_iterator_type, option_iterator_type> p = equal_range(canonicalize(key));
    option_iterator_type it = p.first;
    size_t i = index;
    while (i-- > 0) {
        if (++it == p.second) {
            return "";
        }
    }
    return it->second;
}

bool Configuration::contains(const string& key, const std::size_t n) const {
    pair<option_iterator_type, option_iterator_type> p = equal_range(canonicalize(key));
    option_iterator_type it = p.first;
    size_t i = n;
    while (i-- > 0) {
        if (++it == p.second) {
            return false;
        }
    }
    return true;
}

size_t Configuration::count(const string& key) const {
    pair<option_iterator_type, option_iterator_type> p = equal_range(canonicalize(key));
    option_iterator_type it = p.first;
    size_t i = 0;
    while (++it != p.second) {
        ++i;
    }
    return i;
}


void Configuration::set(const string& key, const string& value) {
    string c_key = canonicalize(key);
    erase(c_key);
    add(c_key, value);
}

void Configuration::erase(const string& key) {
    multimap::erase(canonicalize(key));
}

void Configuration::add(const string& key, const string& value) {
    insert(pair<string, string>(canonicalize(key), value));
}

void Configuration::read_config_line(const string& line) {
    if (line.size() > 2 && line[0] != '#') {
        // Split the line into two pieces on the '=' character
        // The first piece becomes the key, the second becomes the value
        vector<string> kv = split(line, "\n=", 2);
        if (kv.size() < 2) {
            return;
        }
        assert(kv.size() == 2);
        bool replace_config = true;
        string key = kv[0];
        size_t keylen = key.length();
        if (key[keylen-1] == '+') {
            replace_config = false;
            key = trim(key, " \n\t", 0, keylen-1);
        }
        key = canonicalize(key);
        if (replace_config) {
            erase(key);
        }
        // split the value on commas
        vector<string> v = split(kv[1], ",");
        for (vector<string>::iterator it = v.begin(); it != v.end(); it++) {
            string value = trim(*it, " \n\t");
            add(key, value);
        }
    }
}

void Configuration::read_config(istream& in) {
    string line;
    do {
        getline(in, line);
        read_config_line(line);
    } while (!in.eof());
}

Configuration& operator>>(const string& line, Configuration& cc) {
    cc.read_config_line(line);
    return cc;
}

istream& operator>>(istream& in, Configuration& cc) {
    cc.read_config(in);
    return in;
}

ostream& operator<<(ostream& out, const Configuration& cc) {
    string last_key;
    for (multimap<string, string>::const_iterator it = cc.begin(); it != cc.end(); it++) {
        if (last_key == it->first) {
            out << ", " << it->second;
        }
        else {
            if (!last_key.empty()) {
                out << endl;
            }
            out << it->first << " = " << it->second;
        }
        last_key = it->first;
    }
    out << endl;
    return out;
}

ostream& operator<<(ostream& out, Configuration& cc) {
    string last_key;
    for (multimap<string, string>::const_iterator it = cc.begin(); it != cc.end(); it++) {
        if (last_key == it->first) {
            out << ", " << it->second;
        }
        else {
            if (!last_key.empty()) {
                out << endl;
            }
            out << it->first << " = " << it->second;
        }
        last_key = it->first;
    }
    out << endl;
    return out;
}
