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

#pragma once

#ifndef _GSL_EXCEPTION_H_
#define _GSL_EXCEPTION_H_

#include <exception>
#include <gsl/gsl_errno.h>
#include <string>

/**
 * An exception to be thrown when there is an error in the GSL code.
 */
class GSLException : public std::exception {
private:
    std::string _reason;
    std::string _file;
    int _line;
    int _gsl_errno;
    std::string _message;
public:
    GSLException(const char* reason, const char* file, int line, int gsl_errno) throw() :
        _reason(reason), _file(file), _line(line), _gsl_errno(gsl_errno) {
        std::ostringstream s;
        s << "GSL error " << gsl_errno << "(" << gsl_strerror(gsl_errno) << "): " << reason << " at " << file << ":" << line;
        _message = s.str();
    }
    GSLException(const GSLException& e) throw() :
        _reason(e._reason), _file(e._file), _line(e._line), _gsl_errno(e._gsl_errno), _message(e._message) {
    }
    GSLException& operator=(const GSLException& e) throw() {
        _reason = e._reason;
        _file = e._file;
        _line = e._line;
        _gsl_errno = e._gsl_errno;
        _message = e._message;
    }
    ~GSLException() throw() {
    }
    const std::string& reason() const {
        return _reason;
    }
    const std::string& file() const {
        return _file;
    }
    const int line() const {
        return _line;
    }
    const int gsl_errno() const {
        return _gsl_errno;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * GSL error handler function that throws a GSLException.
 */
void gsl_error_throw(const char* reason, const char* file, int line, int gsl_errno) {
    throw GSLException(reason, file, line, gsl_errno);
}

#endif
