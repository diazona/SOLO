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

#include <muParser.h>
#include "hardfactor.h"
#include "../typedefs.h"

class ParsedCompositeHardFactor : public HardFactor {
public:
    ParsedCompositeHardFactor(const string& name, const HardFactor::HardFactorOrder order, const size_t term_count, const HardFactorTerm** terms);
    ParsedCompositeHardFactor(const string& name, const HardFactor::HardFactorOrder order, const HardFactorTermList terms);
    ~ParsedCompositeHardFactor();
    
    const char* get_name() const { return m_name.c_str(); }
    const HardFactorOrder get_order() const { return m_order; }
    const size_t get_term_count() const { return m_term_count; }
    const HardFactorTerm* const* get_terms() const { return m_terms; }
private:
    const bool need_to_free_m_terms;
    const std::string m_name;
    const HardFactorOrder m_order;
    const size_t m_term_count;
    const HardFactorTerm** m_terms;
};

class ParsedHardFactorTerm : public HardFactorTerm {
public:
    ParsedHardFactorTerm(
        const string& name,
        const HardFactorOrder order,
        const IntegrationType* type,
        const string& Fs_real, const string& Fs_imag,
        const string& Fn_real, const string& Fn_imag,
        const string& Fd_real, const string& Fd_imag
    );

    const char* get_name() const { return m_name.c_str(); }
    const IntegrationType* get_type() const { return mp_type; }
    const HardFactorOrder get_order() const { return m_order; }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;

private:
    /**
     * A memory of which IntegrationContext currently provides the variables
     * for Fs_parser.
     * 
     * In order to use a muParser-parsed expression as the hard factor,
     * as this class does, we need to give the Parser a pointer corresponding
     * to each variable that appears in the parsed expression, using
     * Parser::DefineVar. But configuring a Parser with these name-pointer
     * mappings is fairly computationally expensive, and it will slow down the
     * computation a lot if we do it every time Fs() is called.
     * 
     * Fortunately, we don't have to do that. The name-pointer mappings that
     * we need to provide to the Parser are properties of the IntegrationContext
     * begin used. So if we provide the Parser with the pointers from a given
     * instance of IntegrationContext on one call to Fs(), then if the same
     * IntegrationContext instance is passed to Fs() the next time, the same
     * name-pointer mappings are still valid. There's no need to clear out the
     * Parser's variable table and set it up again. Typically, Fs() will be
     * called many times in a row with the same IntegrationContext, so it makes
     * sense to optimize for this case.
     * 
     * Accordingly, at each call to Fs() we save the IntegrationContext provided
     * in this variable. Then at the next call, we can check whether the
     * IntegrationContext passed to Fs() is the same as the one saved. If so,
     * we clear and recreate the name-pointer mappings using the new
     * IntegrationContext, otherwise we do nothing.
     * 
     * This behavior is not really thread-safe, but it will work as long as
     * all calls to Fs() using a given instance of IntegrationContext complete
     * before any calls to Fs() using a different IntegrationContext begin.
     */
    mutable const IntegrationContext* m_ictx_Fs;
    /**
     * A memory of which IntegrationContext currently provides the variables
     * for Fn_parser.
     * 
     * This is necessary for the same reason as m_ictx_Fs.
     */
    mutable const IntegrationContext* m_ictx_Fn;
    /**
     * A memory of which IntegrationContext currently provides the variables
     * for Fd_parser.
     * 
     * This is necessary for the same reason as m_ictx_Fs.
     */
    mutable const IntegrationContext* m_ictx_Fd;
    
    mutable mu::Parser Fs_parser;
    mutable mu::Parser Fn_parser;
    mutable mu::Parser Fd_parser;
    
    const std::string m_name;
    const HardFactorOrder m_order;
    const IntegrationType* mp_type;
};

// This class is NOT thread-safe
class HardFactorParser {
public:
    HardFactorParser();
    void parse_file(const string& filename, bool (*error_handler)(const std::exception& e, const std::string& filename, const size_t line_number) = NULL);
    void parse_line(const string& line);
    HardFactorList get_hard_factors();
    void reset_current_term();

private:
    const bool hard_factor_definition_complete() const;
    const bool hard_factor_definition_empty() const;
    const ParsedHardFactorTerm* create_hard_factor_term();

    HardFactorList hard_factors;
    string Fs_real, Fs_imag, Fn_real, Fn_imag, Fd_real, Fd_imag;
    string name;
    HardFactor::HardFactorOrder order;
    const IntegrationType* type;
    HardFactorRegistry* registry;
};

/**
 * An exception to be thrown when only some of the hard factor properties
 * are defined at the time a ::ParsedHardFactorTerm is created.
 * 
 * This exception gets thrown when you have something like this
 * in a hard factor definition file:
 * 
 *     name = foo
 *     type = none
 *     Fs real = x
 *     type = momentum1
 * 
 * The parser reads a name, a type, and a real Fs, then it reads another name
 * and tries to start a new hard factor specification with that line. But there
 * are no definitions for Fs imag, Fn real, Fn imag, Fd real, Fd imag, or
 * order. So there's not enough information to make a hard factor term.
 */
class IncompleteHardFactorDefinitionException : public exception {
private:
    string _message;
public:
    /**
     * Constructs an instance of the exception.
     *
     * @param message a descriptive human-readable message
     */
    IncompleteHardFactorDefinitionException(const string& message = "") throw() {
        if (message.empty()) {
            _message = "Incomplete hard factor definition";
        }
        else {
            _message = "Incomplete hard factor definition: " + message;
        }
    }
    ~IncompleteHardFactorDefinitionException() throw() {};
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * An exception to be thrown when a hard factor property has an invalid value.
 */
class InvalidHardFactorDefinitionException : public std::exception {
private:
    std::string _message;
public:
    const std::string line;
    const std::string key;
    const std::string value;
    /**
     * Constructs an instance of the exception.
     *
     * @param message a descriptive human-readable message
     */
    InvalidHardFactorDefinitionException(const std::string& line, const std::string& key, const std::string& value, const std::string& message = "") throw() {
        std::ostringstream oss;
        if (message.empty()) {
            oss << "Invalid hard factor definition: ";
        }
        else {
            oss << message << " ";
        }
        oss << key << " = " << value;
        _message = oss.str();
    }
    virtual ~InvalidHardFactorDefinitionException() throw() {};
    const char* what() const throw() {
        return _message.c_str();
    }
};

