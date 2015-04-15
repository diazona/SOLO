/*
 * Part of SOLO
 *
 * Copyright 2014 David Zaslavsky
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

#pragma once

#include <stdexcept>
#include <list>
#include <map>
#include <string>

/**
 * A template class used by ::HardFactorRegistry to store hard
 * factors and groups.
 */
template<typename T>
class category_map : public std::map<const std::string, std::list<std::pair<const std::string, T> > > {
private:
    typedef typename std::map<const std::string, std::list<std::pair<const std::string, T> > > super;
public:
    typedef std::pair<const std::string, T> pair_type;
    typedef std::list<pair_type> list_type;

    category_map() : super() {}
    category_map(const category_map<T>& other) : super(other) {}
    category_map<T>& operator=(const category_map<T>& other) {
        super::operator=(other);
        return *this;
    }
    ~category_map() {}

    static void normalize(std::string& s) {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    }

    void add(const std::string& name, const std::string& implementation, const T t) {
        std::string l_name(name);
        std::string l_impl(implementation);
        normalize(l_name);
        normalize(l_impl);

        list_type& the_list = super::operator[](l_name);
        for (typename list_type::iterator it = the_list.begin(); it != the_list.end(); it++) {
            if (it->first == l_impl) {
                it->second = t;
                return;
            }
        }
        pair_type the_pair(l_impl, t);
        the_list.push_back(the_pair);
    }

    const T& get(const std::string& name, const std::string& implementation) const {
        std::string l_name(name);
        std::string l_impl(implementation);
        normalize(l_name);
        normalize(l_impl);

        typename super::const_iterator map_iterator = super::find(l_name);
        if (map_iterator != super::end()) {
            const list_type& the_list = map_iterator->second;
            for (typename list_type::const_iterator it = the_list.begin(); it != the_list.end(); it++) {
                if (it->first == l_impl) {
                    return it->second;
                }
            }
        }
        throw std::out_of_range("no element with name " + name + " and implementation " + implementation);
    }

    const T& get(const std::string& name, const bool first = true) const {
        std::string l_name(name);
        normalize(l_name);

        typename super::const_iterator map_iterator = super::find(l_name);
        if (map_iterator != super::end()) {
            const list_type& the_list = map_iterator->second;
            if (!the_list.empty()) {
                const pair_type& pr = first ? the_list.front() : the_list.back();
                const T& p = pr.second;
                return p;
            }
        }
        throw std::out_of_range("no element with name " + name);
    }
};
