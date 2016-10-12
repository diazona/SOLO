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

#include <list>
#include <string>
#include <muParser.h>
#include "hardfactor.h"

/**
 * A ::HardFactorTerm subclass which represents formulas parsed from text.
 */
class ParsedHardFactorTerm : public HardFactorTerm {
public:
    ParsedHardFactorTerm(
        const std::string& name,
        const std::string& implementation,
        const HardFactor::HardFactorOrder order,
        const IntegrationRegion* region,
        const Modifiers& modifiers,
        const std::string& Fs_real, const std::string& Fs_imag,
        const std::string& Fn_real, const std::string& Fn_imag,
        const std::string& Fd_real, const std::string& Fd_imag,
        const std::list<std::pair<std::string, std::string> >& variable_list
    );
    ParsedHardFactorTerm(
        const std::string& name,
        const std::string& implementation,
        const HardFactor::HardFactorOrder order,
        const CoreIntegrationRegion& core,
        const std::vector<const AuxiliaryIntegrationRegion*> subregions,
        const Modifiers& modifiers,
        const std::string& Fs_real, const std::string& Fs_imag,
        const std::string& Fn_real, const std::string& Fn_imag,
        const std::string& Fd_real, const std::string& Fd_imag,
        const std::list<std::pair<std::string, std::string> >& variable_list
    );
    ~ParsedHardFactorTerm();

    const char* get_name() const { return m_name.c_str(); }
    const char* get_implementation() const { return m_implementation.c_str(); }
    const IntegrationRegion* get_integration() const { return mp_region; }
    HardFactorOrder get_order() const { return m_order; }
    const Modifiers& get_modifiers() const { return m_modifiers; }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;

    const std::string Fs_expr() const;
    const std::string Fn_expr() const;
    const std::string Fd_expr() const;

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
     * IntegrationContext passed to Fs() is the same as the one saved. If not,
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
    const std::string m_implementation;
    const HardFactorOrder m_order;
    const Modifiers m_modifiers;

    const IntegrationRegion* mp_region;
    const bool m_free_region;

    size_t aux_variable_count;
    double* aux_variable_storage;

    void init_term(
        const std::string& Fs_real, const std::string& Fs_imag,
        const std::string& Fn_real, const std::string& Fn_imag,
        const std::string& Fd_real, const std::string& Fd_imag,
        const std::list<std::pair<std::string, std::string> >& variable_list
    );
    void init_parser(mu::Parser& parser, mu::varmap_type& all_used_variables, const string& aux_variables, const string& real_expr, const string& imag_expr, const char* debug_message);
#ifndef NDEBUG
    void print_parser_info(const char* message, mu::Parser& parser, const double* real, const double* imag) const;
    void print_parser_info(const char* message, mu::Parser& parser) const;
# endif
};

/**
 * A ::HardFactor representing a sum of multiple terms, where the specification
 * of which terms has been parsed from text.
 */
class ParsedCompositeHardFactor : public HardFactor {
public:
    ParsedCompositeHardFactor(const string& name, const string& implementation, const HardFactor::HardFactorOrder order, const size_t term_count, const HardFactorTerm** terms);
    ParsedCompositeHardFactor(const string& name, const string& implementation, const HardFactor::HardFactorOrder order, const HardFactorTermList terms);
    ~ParsedCompositeHardFactor();

    const char* get_name() const { return m_name.c_str(); }
    const char* get_implementation() const { return m_implementation.c_str(); }
    HardFactorOrder get_order() const { return m_order; }
    size_t get_term_count() const { return m_term_count; }
    const HardFactorTerm* const* get_terms() const { return m_terms; }
private:
    const bool need_to_free_m_terms;
    const std::string m_name;
    const std::string m_implementation;
    const HardFactorOrder m_order;
    const size_t m_term_count;
    const HardFactorTerm** m_terms;
};

/**
 * The object that parses hard factors from text specifications.
 */
class HardFactorParser {
public:
    HardFactorParser(HardFactorRegistry& registry);
    ~HardFactorParser();

    /**
     * Parses a file containing hard factor specifications. The parser recognizes
     * three kinds of content:
     *
     * - Elementary hard factor specifications, which contain a name, an integration
     *   region specification, optional evaluation modifiers, as well as formulas
     *   to be used to evaluate the parts of the hard factor. Each of these is given
     *   on a line in `key = value` format. For example:
     *
     *     name = testhf1
     *     integration = nlo clipped * cartesian position
     *     modifiers = exact, divide xi
     *     Fs real = 1 + xi2
     *     Fs imag = 1 - xi2
     *     Fn real = 1/z
     *     Fn imag = 1/z
     *     Fd real = r2/(s2*t2)
     *     Fd imag = r2/(s2*t2)
     *
     *   The variables that can be used in the expressions are those stored in
     *   ::IntegrationContext and ::Context.
     *
     *   Technically speaking, these correspond to instances of ::HardFactorTerm,
     *   but in most cases they are only used as instances of ::HardFactor.
     * - Composite hard factor specifications, which contain a name, an optional
     *   implementaion code, and a list of elementary hard factor specifications
     *   which constitute the composite hard factor. The composite hard factor
     *   represents the sum of the elementary hard factors that it contains. This
     *   specification is given in the format
     *
     *     hfname.hfimpl = hfname.hfimpl,hfname.hfimpl
     *
     *   with the name and implementation of the composite hard factor occuring
     *   as the key (before the equal sign), and the names and implementations of
     *   the elementary hard factors that make it up as the value, separated by
     *   commas. Implementation codes are optional in all cases (though bear in mind
     *   that things can get weird if multiple hard factors are defined with the same
     *   name but one of them doesn't have an implementation code). Examples:
     *
     *     testchf.foo = testhf1.m,testhf2.m,testhf3
     *     testchf.bar = testhf1.r,testhf2.r,testhf3
     *
     *   Once a composite hard factor is defined in this way, it can be used just like
     *   an elementary hard factor (i.e. as an instance of ::HardFactor).
     * - Hard factor group specifications, which contain a label and a list of
     *   (elementary or composite) hard factor specifications. These are given as
     *
     *     label.hfgimpl:hfname.hfimpl,hfname.hfimpl,hfname.hfimpl
     *
     *   Groups differ from composite hard factors in that their constitutent hard factors
     *   can be separated by the program and calculated individually. Groups are not
     *   actually parsed immediately; instead, a line that constitutes a group
     *   specification is stored until ::flush_groups() is called, and that does
     *   the actual parsing.
     *
     * Each time a hard factor definition (elementary or composite) is read from
     * the file, this calls the callback registered with ::set_hard_factor_callback,
     * and each time a hard factor group definition is read from the file, it
     * calls the callback registered with ::set_hard_factor_group_callback.
     *
     * This method basically just executes ::parse_line for each line of the file
     * in order, with a bit of dressing which calls the error handler registered
     * using ::set_error_handler in case `parse_line` ever throws an error.
     *
     * @param[in] filename the name of the file to parse
     */
    void parse_file(const std::string& filename);

    /**
     * Parses a single line of a file. This is the method that actually does the
     * parsing, but most client code never calls it directly. To parse a file you
     * should usually use ::parse_file instead. This method is only exposed as
     * part of the public interface in case some client code wants to do something
     * weird that `parse_file` can't handle (for example, parsing a stream).
     *
     * @param[in] line the line
     */
    void parse_line(const std::string& line);

    /**
     * Attempts to parse a string as a hard factor group specification.
     * The specification passed should take the form
     *
     *     name:hfname.hfimpl,hfname.hfimpl,hfname.hfimpl,...
     *
     * where
     *
     * - `name` is the descriptive name of the hard factor group,  which can
     *   be used to refer to it on the command line of `oneloopcalc`. The name
     *   is optional; if one is not provided, the hard factor group's name will
     *   simply be the full text of the specification.
     * - `hfname.hfimpl` is a hard factor specification, as described in the
     *   documentation for ::parse_hardfactor.
     *
     * This method is used internally by ::flush_groups if a hard factor group
     * specification is encountered while parsing, or it can also be called
     * by external code. It always returns a newly constructed ::HardFactorGroup
     * object. The new object is _not_ added to the ::HardFactorRegistry by this
     * method, though when this is called from `flush_groups()`, the returned pointer
     * will be added to the registry.
     *
     * @return the group object, or `NULL` if the parsing failed to produce
     * a valid group
     */
    const HardFactorGroup* parse_hard_factor_group(const std::string& spec);

    /**
     * Splits a string representing a hard factor specification into a name
     * and implementation. The specification passed should take the form
     *
     *     name.implementation
     *
     * or just
     *
     *     name
     *
     * where
     *
     * - `name` is the descriptive name of the hard factor.
     * - `implementation` identifies which implementation of the hard factor
     *   is desired. Sometimes there are multiple expressions that can be
     *   used to calculate the same thing; the implementation code is a way
     *   to distinguish between them. For example, `m` or `momentum` means
     *   to use a momentum-space expression, `r` or `radial` for a radial
     *   position space expression, and `p` or `position` for a Cartesian
     *   position space expression. The implementation code is optional.
     *
     * This function makes no attempt to verify that either `name` or
     * `implementation` actually refers to a real hard factor. To do that,
     * use ::parse_hardfactor.
     *
     * @param[in] spec the hard factor specification
     * @param[out] name the name part of the specification
     * @param[out] implementation the implementation part of the specification,
     * or the empty string if no implementation part is provided
     */
    static void split_hardfactor(const std::string& spec, std::string& name, std::string& implementation);

    /**
     * This method goes over all stored group specifications and parses them
     * using whatever hard factors are currently entered into the registry.
     *
     * If a HardFactorParser is deleted or goes out of scope,
     * `flush_groups()` will be automatically called from the destructor.
     */
    void flush_groups();

    void set_hard_factor_callback(void (*callback)(const HardFactor& hf));
    void set_hard_factor_group_callback(void (*callback)(const HardFactorGroup& hfg));
    void set_error_handler(bool (*error_handler)(const std::exception& e, const std::string& filename, const size_t line_number));

    HardFactorRegistry& registry;
private:
    // these shouldn't be called
    HardFactorParser(HardFactorParser& other) : registry(other.registry) {}
    HardFactorParser& operator=(HardFactorParser& other) { registry = other.registry; return *this; }

    bool hard_factor_definition_complete() const;
    bool hard_factor_definition_empty() const;
    const ParsedHardFactorTerm* create_hard_factor_term();
    void reset_current_term();
    bool unflushed_groups();

    const ParsedCompositeHardFactor* parse_composite_hard_factor(const string& key, const string& value);
    void parse_variable_definition(const string& key, const string& value);

    HardFactorList hard_factors;
    std::list<std::string> unparsed_hard_factor_group_specs;
    std::list<std::pair<std::string, std::string> > variable_definitions;
    std::string Fs_real, Fs_imag, Fn_real, Fn_imag, Fd_real, Fd_imag;
    std::string name, implementation, default_implementation;
    HardFactor::HardFactorOrder order;
    Modifiers modifiers;
    CoreIntegrationRegion* core;
    vector<const AuxiliaryIntegrationRegion*> subregions;

    // callbacks
    void (*hard_factor_callback)(const HardFactor& hf);
    void (*hard_factor_group_callback)(const HardFactorGroup& hfg);
    bool (*error_handler)(const std::exception& e, const std::string& filename, const size_t line_number);
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
    IncompleteHardFactorDefinitionException(const string& message = "") throw();
    ~IncompleteHardFactorDefinitionException() throw() {};
    const char* what() const throw();
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
    InvalidHardFactorDefinitionException(const std::string& line, const std::string& key, const std::string& value, const std::string& message = "") throw();
    virtual ~InvalidHardFactorDefinitionException() throw() {};
    const char* what() const throw();
};

/**
 * An exception to be thrown when a hard factor specification fails to be
 * parsed for some reason. The reason is indicated by the message.
 */
class InvalidHardFactorSpecException : public std::exception {
private:
    std::string _message;
public:
    std::string hfspec;
    /**
     * Constructs an instance of the exception.
     *
     * @param hfspec the hard factor specification or part of a specification
     *  that caused the error
     * @param message a descriptive human-readable message indicating why hfspec
     *  could not be parsed
     */
    InvalidHardFactorSpecException(const std::string& hfspec, const std::string& message) throw();
    InvalidHardFactorSpecException(const InvalidHardFactorSpecException& other) throw() : _message(other._message) {}
    ~InvalidHardFactorSpecException() throw() {}
    void operator=(const InvalidHardFactorSpecException& other);
    const char* what() const throw();
};
