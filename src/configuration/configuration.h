/*
 * Part of SOLO
 *
 * Copyright 2015 David Zaslavsky
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

/**
 * Configuration file parsers are a dime a dozen, but I don't know of one that
 * meets the requirements I have for SOLO. Specifically, it needs to allow
 * keys to contain a list of values which can be added to or replaced, not
 * just a single value. It also needs to ignore anything that doesn't look like
 * a configuration parameter. However, I don't need sections or structures,
 * and I don't need type safety in configuration options. Just a dumb
 * (string -> list of string) container with some parsing capability.
 *
 * So I'm writing it myself.
 */

#pragma once

#include <map>
#include <string>

/**
 * A repository for settings. A Configuration object is able to read a
 * configuration file in
 *  key = value
 * format, and store all the settings read.
 */
class Configuration : public std::multimap<std::string, std::string> {
public:
    /**
     * Construct an empty Configuration.
     */
    Configuration(std::string (*canonicalize)(const std::string&) = NULL);

    /**
     * Construct a Configuration and initialize it with settings
     * read from the named file.
     */
    Configuration(const std::string& filename, std::string (*canonicalize)(const std::string&) = NULL);

    /**
     * Removes all settings with the given key.
     */
    void erase(const std::string& key);
    /**
     * Returns the value of a setting, if any, or the empty string if unset.
     */
    std::string get(const std::string& key, const size_t index = 0) const;
    /**
     * Returns whether the given key has a value (or at least n values)
     * associated with it.
     */
    bool contains(const std::string& key, const size_t n = 1) const;
    /**
     * Returns the number of values associated with the given key.
     */
    size_t count(const std::string& key) const;
    /**
     * Add a setting, replacing any existing settings with the same key.
     */
    void set(const std::string& key, const std::string& value);
    /**
     * Add a setting. Any existing settings with the same key are left alone;
     * in this case there will be multiple values with that key.
     */
    void add(const std::string& key, const std::string& value);

    /**
     * Process a string representing one line of a config file (i.e. one setting)
     */
    void read_config_line(const std::string& line);

    /**
     * Read a config file, or something in an equivalent format, from an input stream
     * and add the settings to the Configuration.
     */
    void read_config(std::istream& in);

    friend Configuration& operator>>(std::string& line, Configuration& cc);
    friend std::istream& operator>>(std::istream& in, Configuration& cc);
    friend std::ostream& operator<<(std::ostream& out, Configuration& cc);

private:
    std::string (*canonicalize)(const std::string&);
};

/**
 * Allows adding a setting to a Configuration using >> notation.
 */
Configuration& operator>>(std::string& line, Configuration& cc);
/**
 * Allows reading a Configuration in from a stream using >> notation.
 */
std::istream& operator>>(std::istream& in, Configuration& cc);
/**
 * Allows writing a Configuration out to a stream using << notation.
 *
 * What is written out is just the list of key-value pairs. The output
 * could be read in to reconstruct the Configuration.
 */
std::ostream& operator<<(std::ostream& out, Configuration& cc);
/**
 * Allows writing a Configuration out to a stream using << notation.
 *
 * What is written out is just the list of key-value pairs. The output
 * could be read in to reconstruct the Configuration.
 */
std::ostream& operator<<(std::ostream& out, const Configuration& cc);
