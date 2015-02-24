/*
 * Part of oneloopcalc
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

#include <cassert>
#include <fstream>
#include <istream>
#include <sstream>
#include <muParser.h>
#include <gsl/gsl_math.h>
#include "gsl_mu.h"
#include "hardfactor.h"
#include "hardfactor_parsed.h"
#include "hardfactor_parser.h"
#include "../typedefs.h"
#include "../utils/utils.h"

using mu::Parser;
using mu::value_type;
using mu::valmap_type;
using mu::varmap_type;

ParsedCompositeHardFactor::ParsedCompositeHardFactor(const string& name, const HardFactor::HardFactorOrder order, const size_t term_count, const HardFactorTerm** terms) :
  HardFactor(), 
  m_name(name), m_order(order), m_term_count(term_count), m_terms(terms), need_to_free_m_terms(false) {
}

ParsedCompositeHardFactor::ParsedCompositeHardFactor(const string& name, const HardFactor::HardFactorOrder order, const HardFactorTermList terms) :
  HardFactor(), 
  m_name(name), m_order(order), m_term_count(terms.size()), m_terms(new const HardFactorTerm*[m_term_count]), need_to_free_m_terms(true) {
    const HardFactorTerm** m_terms_end = copy(terms.begin(), terms.end(), m_terms);
    assert(m_terms_end - m_terms == m_term_count);
}

ParsedCompositeHardFactor::~ParsedCompositeHardFactor() {
    if (need_to_free_m_terms) {
        delete[] m_terms;
    }
}


/**
 * Return the string, or "0" if the string is empty
 */
static const string e0(const string& s) {
    return s.empty() ? "0" : s;
}

#ifndef NDEBUG
// the only way to set this to 1 is through a debugger
static volatile int deep_debugging = 0;

void ParsedHardFactorTerm::print_parser_info(const char* message, Parser& parser, const double* real, const double* imag) const {
    if (deep_debugging > 0) {
        cerr << "(" << this << ")" << get_name();
        if (message != NULL) {
            cerr << " " << message;
        }
        cerr << ": " << parser.GetExpr();
        if (real != NULL && imag != NULL) {
            cerr << " = " << *real << " + " << *imag << "i";
        }
        cerr << endl;
    }
}

void ParsedHardFactorTerm::print_parser_info(const char* message, Parser& parser) const {
    print_parser_info(message, parser, NULL, NULL);
}
# endif

// until I get a better solution, whether altering muParser
// or some sort of machine code voodoo
static GluonDistribution* gdist = NULL;

value_type gluon_distribution_F(const value_type k2, const value_type Y) {
    assert(gdist != NULL);
    return gdist->F(k2, Y);
}
value_type gluon_distribution_S2(const value_type r2, const value_type Y) {
    assert(gdist != NULL);
    return gdist->S2(r2, Y);
}
value_type gluon_distribution_S4(const value_type r2, const value_type s2, const value_type t2, const value_type Y) {
    assert(gdist != NULL);
    return gdist->S4(r2, s2, t2, Y);
}

/**
 * Dot product of 2D vectors.
 */
value_type dot2(const value_type a1, const value_type a2, const value_type b1, const value_type b2) {
    return a1*b1 + a2*b2;
}
/**
 * Square of a 2D vector.
 */
value_type square2(const value_type a1, const value_type a2) {
    return dot2(a1, a2, a1, a2);
}
/**
 * Norm (magnitude) of a 2D vector.
 */
value_type norm2(const value_type a1, const value_type a2) {
    return gsl_hypot(a1, a2);
}

void ParsedHardFactorTerm::init_parser(Parser& parser, varmap_type& all_variables, const string& real_expr, const string& imag_expr, const char* debug_message) {
    parser.SetExpr(e0(real_expr) + "," + e0(imag_expr));
    mu_load_gsl(parser);
    parser.DefineFun("F", gluon_distribution_F);
    parser.DefineFun("S2", gluon_distribution_S2);
    parser.DefineFun("S4", gluon_distribution_S4);
    parser.DefineFun("dot", dot2);
    parser.DefineFun("square", square2);
    parser.DefineFun("norm", norm2);
#ifndef NDEBUG
    print_parser_info(debug_message, parser);
# endif
    // use GetUsedVar() to trigger an attempt to parse the expression
    // we want any parsing errors to be revealed now, not when evaluating expressions
    varmap_type variables_from_parser = parser.GetUsedVar();
    all_variables.insert(variables_from_parser.begin(), variables_from_parser.end());
}

ParsedHardFactorTerm::ParsedHardFactorTerm(
    const string& name,
    const HardFactor::HardFactorOrder order,
    const IntegrationType* type,
    const string& Fs_real, const string& Fs_imag,
    const string& Fn_real, const string& Fn_imag,
    const string& Fd_real, const string& Fd_imag) :
  m_name(name),
  m_order(order),
  mp_type(type) {
    varmap_type all_variables;
    init_parser(Fs_parser, all_variables, Fs_real, Fs_imag, "Fs");
    init_parser(Fn_parser, all_variables, Fn_real, Fn_imag, "Fn");
    init_parser(Fd_parser, all_variables, Fd_real, Fd_imag, "Fd");

    // make sure only existing variables are used
#define process(var) all_variables.erase(#var);
#include "../integration/ictx_var_list.inc"
#include "../configuration/ctx_var_list.inc"
#undef process
    if (!all_variables.empty()) {
        ostringstream oss;
        oss << "Invalid variables used:";
        for (varmap_type::const_iterator it = all_variables.begin(); it != all_variables.end(); it++) {
            oss << " " << it->first;
        }
        throw mu::ParserError(oss.str());
    }
}

void define_variables(Parser& parser, const IntegrationContext* ictx) {
    /* An alternative implementation would be to have an IntegrationContext
     * instance in ParsedHardFactorTerm itself - not just a pointer, a
     * concrete instance - and then this function would just copy the
     * passed IntegrationContext to the member IntegrationContext.
     * @code m_ictx = ictx;
     * But then, if the passed IntegrationContext has its values changed
     * elsewhere in the program, those changes won't be reflected in the
     * Parser. That's why I don't do that.
     */
#define process(var) parser.DefineVar(#var, const_cast<value_type*>(&(ictx->var)));
#include "../integration/ictx_var_list.inc"
#undef process
#define process(var) parser.DefineVar(#var, const_cast<value_type*>(&(ictx->ctx->var)));
#include "../configuration/ctx_var_list.inc"
#undef process
    // and now some aliases
    parser.DefineVar("A", const_cast<value_type*>(&(ictx->ctx->mass_number)));
    parser.DefineVar("c", const_cast<value_type*>(&(ictx->ctx->centrality)));
    
    /* The last "variables" we need on a per-IntegrationContext basis are
     * the gluon distribution functions. These functions need access to
     * ictx->gdist so they can call ictx->gdist->F(k2, Y) or the like.
     * However, when I add a two-argument function to the Parser lexicon,
     * muParser's API doesn't let me pass &(ictx->gdist) to the underlying
     * C function. The simplest way to work around this is to use a global
     * variable to store the gluon distribution, which is fine as long as
     * there is only ever one gluon distribution used. If the program ever
     * changes to use more than one gluon distribution per run, then this
     * static global variable will be insufficient and I'll need to come up
     * with some more sophisticated plan.
     */
    if (gdist == NULL) {
        gdist = ictx->ctx->gdist;
    }
    else {
        /* I think it's safe to make this an assertion (which gets stripped
         * in a release build) rather than a simple if statement (which would
         * still be included in a release build) because in the current version
         * of the program there's no way this assertion should ever be false.
         * The program doesn't allow you to define multiple gluon distributions.
         */
        assert(gdist == ictx->ctx->gdist);
    }
    
    // just in case ictx->ctx->gdist is NULL or something,  I dunno
    assert(gdist != NULL);
}

void evaluate_hard_factor(Parser& parser, double* real, double* imag) {
    int number_of_values;
    value_type* values;
    values = parser.Eval(number_of_values);
    if (number_of_values != 2) {
        throw mu::ParserError("invalid number of values");
        return;
    }
    *real = values[0];
    *imag = values[1];
}

void ParsedHardFactorTerm::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx != m_ictx_Fs) {
        m_ictx_Fs = ictx;
        define_variables(Fs_parser, ictx);
    }
    evaluate_hard_factor(Fs_parser, real, imag);
#ifndef NDEBUG
    print_parser_info("Fs", Fs_parser, real, imag);
# endif
}

void ParsedHardFactorTerm::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx != m_ictx_Fn) {
        m_ictx_Fn = ictx;
        define_variables(Fn_parser, ictx);
    }
    evaluate_hard_factor(Fn_parser, real, imag);
#ifndef NDEBUG
    print_parser_info("Fn", Fn_parser, real, imag);
# endif
}

void ParsedHardFactorTerm::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx != m_ictx_Fd) {
        m_ictx_Fd = ictx;
        define_variables(Fd_parser, ictx);
    }
    evaluate_hard_factor(Fd_parser, real, imag);
#ifndef NDEBUG
    print_parser_info("Fd", Fd_parser, real, imag);
# endif
}

const string ParsedHardFactorTerm::Fs_expr() const {
    return Fs_parser.GetExpr();
}

const string ParsedHardFactorTerm::Fn_expr() const {
    return Fn_parser.GetExpr();
}

const string ParsedHardFactorTerm::Fd_expr() const {
    return Fd_parser.GetExpr();
}

using std::ifstream;
using std::vector;

static HardFactor::HardFactorOrder sentinel = static_cast<HardFactor::HardFactorOrder>(-1);

HardFactorParser::HardFactorParser() : order(sentinel), type(NULL), registry(NULL) {
}

HardFactorList HardFactorParser::get_hard_factors() {
    return hard_factors;
}


void HardFactorParser::parse_line(const string& line) {
    if (line.length() == 0) {
        return;
    }
    /* A hard factor definition consists of a name, an order
        * (LO or NLO), an integration type, and six functions
        * giving the real and imaginary components of Fs, Fn, Fd.
        * Each part of the definition is given as name = value.
        */
    vector<string> parts = split(line, "=", 2);
    assert(parts.size() > 0);
    if (parts.size() != 2) {
        throw InvalidHardFactorDefinitionException(line, parts[0], parts.size() > 1 ? parts[1] : "", "Incomplete or malformed property");
    }
    parts[0] = trim(parts[0]);
    parts[1] = trim(parts[1]);
    if (parts[0] == "name") {
        if (!name.empty()) {
            create_hard_factor_term();
        }
        name = parts[1];
    }
    else if (parts[0] == "order") {
        if (order != sentinel) {
            create_hard_factor_term();
        }
        if (parts[1] == "lo") {
            order = HardFactor::LO;
        }
        else if (parts[1] == "nlo") {
            order = HardFactor::NLO;
        }
        /* no support for mixed-order hard factors because
            * this function only creates HardFactorTerms which
            * can only be LO or NLO */
        else if (parts[1] == "mixed") {
            throw InvalidHardFactorDefinitionException(line, parts[0], parts[1], "Mixed-order hard factors not supported:");
        }
        else {
            throw InvalidHardFactorDefinitionException(line, parts[0], parts[1]);
        }
    }
    else if (parts[0] == "type") {
        if (type != NULL) {
            create_hard_factor_term();
        }
        if (parts[1] == "none") {
            type = &momentum::none;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "momentum1") {
            type = &momentum::momentum1;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "momentum2") {
            type = &momentum::momentum2;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "momentum3") {
            type = &momentum::momentum3;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "radial momentum1") {
            type = &momentum::radialmomentum1;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "radial momentum2") {
            type = &momentum::radialmomentum2;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "radial momentum3") {
            type = &momentum::radialmomentum3;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "momentumxip1") {
            type = &momentum::momentumxip1;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "momentumxip2") {
            type = &momentum::momentumxip2;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "qlim") {
            type = &momentum::qlim;
            registry = momentum::registry::get_instance();
        }
        else if (parts[1] == "cartesian dipole") {
            type = &position::dipole;
            registry = position::registry::get_instance();
        }
        else if (parts[1] == "cartesian quadrupole") {
            type = &position::quadrupole;
            registry = position::registry::get_instance();
        }
        else if (parts[1] == "radial dipole") {
            type = &radial::dipole;
            registry = radial::registry::get_instance();
        }
        else if (parts[1] == "radial quadrupole") {
            type = &radial::quadrupole;
            registry = radial::registry::get_instance();
        }
        else if (parts[1] == "radial rescaled dipole") {
            type = &radial::rescaled_dipole;
            registry = radial::registry::get_instance();
        }
        else if (parts[1] == "radial rescaled quadrupole") {
            type = &radial::rescaled_quadrupole;
            registry = radial::registry::get_instance();
        }
        else {
            throw InvalidHardFactorDefinitionException(line, parts[0], parts[1], "Unknown integration type:");
        }
    }
    else if (parts[0] == "Fs real") {
        if (!Fs_real.empty()) {
            create_hard_factor_term();
        }
        Fs_real = parts[1];
    }
    else if (parts[0] == "Fs imag") {
        if (!Fs_imag.empty()) {
            create_hard_factor_term();
        }
        Fs_imag = parts[1];
    }
    else if (parts[0] == "Fn real") {
        if (!Fn_real.empty()) {
            create_hard_factor_term();
        }
        Fn_real = parts[1];
    }
    else if (parts[0] == "Fn imag") {
        if (!Fn_imag.empty()) {
            create_hard_factor_term();
        }
        Fn_imag = parts[1];
    }
    else if (parts[0] == "Fd real") {
        if (!Fd_real.empty()) {
            create_hard_factor_term();
        }
        Fd_real = parts[1];
    }
    else if (parts[0] == "Fd imag") {
        if (!Fd_imag.empty()) {
            create_hard_factor_term();
        }
        Fd_imag = parts[1];
    }
    /* Anything other than those keys, if it doesn't contain spaces, is
     * taken to represent a specification of a hard factor which
     * contains multiple terms.
     */
    else if (parts[0].find_first_of(" \t\r\n") == string::npos) {
        create_hard_factor_term();
        vector<string> term_labels = split(parts[1], ",");
        if (term_labels.empty()) {
            throw InvalidHardFactorSpecException(parts[0], "no hard factor terms provided in definition");
        }
        HardFactorTermList hftlist;
        HardFactor::HardFactorOrder order = sentinel;
        for (vector<string>::const_iterator it = term_labels.begin(); it != term_labels.end(); it++) {
            string s = trim(*it);
            if (s.empty()) {
                continue;
            }
            // possible future enhancement: handle the case where s names
            // a composite hard factor by extracting the terms and adding them
            // individually
            const HardFactor* hf = parse_hardfactor(s);
            if (hf == NULL) {
                throw InvalidHardFactorSpecException(s, "no hard factor with that name");
            }
            const HardFactorTerm* hft = dynamic_cast<const HardFactorTerm*>(hf);
            if (hft == NULL) {
                throw InvalidHardFactorSpecException(s, "can't handle composite hard factors");
            }
            hftlist.push_back(hft);
            if (order == sentinel) {
                order = hft->get_order();
            }
            else if (order != HardFactor::MIXED && order != hft->get_order()) {
                order = HardFactor::MIXED;
            }
        }
        if (hftlist.empty()) {
            throw InvalidHardFactorSpecException(parts[1], "no valid hard factors provided in definition");
        }
        string name;
        HardFactorRegistry* registry = parse_one_hardfactor_spec(parts[0], name, false);
        if (registry == NULL) {
            throw InvalidHardFactorSpecException(parts[0], "couldn't identify registry from hard factor name");
        }
        ParsedCompositeHardFactor* hf = new ParsedCompositeHardFactor(name, order, hftlist);
        registry->add_hard_factor(hf, true);
        hard_factors.push_back(hf);
    }
    else {
        throw InvalidHardFactorDefinitionException(line, parts[0], parts[1], "Unknown property:");
    }
}


void HardFactorParser::parse_file(const string& filename, bool (*error_handler)(const exception&, const string&, const size_t line_number)) {
    ifstream hfinput(filename.c_str());
    string line;
    size_t i;

    for (i = 0, getline(hfinput, line); hfinput; i++, getline(hfinput, line)) {
        try {
            parse_line(line);
        }
        catch (const InvalidHardFactorDefinitionException& e) {
            if (error_handler == NULL || error_handler(e, line, i)) {
                throw;
            }
        }
        catch (const IncompleteHardFactorDefinitionException& e) {
            if (error_handler == NULL || error_handler(e, line, i)) {
                throw;
            }
        }
    }
    if (!hard_factor_definition_empty()) {
        create_hard_factor_term();
    }
}

const bool HardFactorParser::hard_factor_definition_empty() const {
    return
      type == NULL &&
      registry == NULL &&
      order == sentinel &&
      name.empty() &&
       (   Fs_real.empty() && Fs_imag.empty()
        && Fn_real.empty() && Fn_imag.empty()
        && Fd_real.empty() && Fd_imag.empty());
}

const bool HardFactorParser::hard_factor_definition_complete() const {
    return
      type != NULL &&
      registry != NULL &&
      order != sentinel &&
      !name.empty() &&
      !(   Fs_real.empty() && Fs_imag.empty()
        && Fn_real.empty() && Fn_imag.empty()
        && Fd_real.empty() && Fd_imag.empty());
}

const ParsedHardFactorTerm* HardFactorParser::create_hard_factor_term() {
    if (!hard_factor_definition_complete()) {
        throw IncompleteHardFactorDefinitionException();
    }
    ParsedHardFactorTerm* hf = new ParsedHardFactorTerm(name, order, type, Fs_real, Fs_imag, Fn_real, Fn_imag, Fd_real, Fd_imag);
    registry->add_hard_factor(hf, true);
    hard_factors.push_back(hf);
    reset_current_term();
    return hf;
}

void HardFactorParser::reset_current_term() {
    // We normally shouldn't be resetting when there is an incomplete definition
    assert(hard_factor_definition_complete());
    name.erase();
    order = sentinel;
    type = NULL;
    registry = NULL;
    Fs_real.erase();
    Fs_imag.erase();
    Fn_real.erase();
    Fn_imag.erase();
    Fd_real.erase();
    Fd_imag.erase();
    assert(hard_factor_definition_empty());
}

IncompleteHardFactorDefinitionException::IncompleteHardFactorDefinitionException(const string& message) throw() {
    if (message.empty()) {
        _message = "Incomplete hard factor definition";
    }
    else {
        _message = "Incomplete hard factor definition: " + message;
    }
}

const char* IncompleteHardFactorDefinitionException::what() const throw() {
    return _message.c_str();
}

InvalidHardFactorDefinitionException::InvalidHardFactorDefinitionException(const string& line, const string& key, const string& value, const string& message) throw() {
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

const char* InvalidHardFactorDefinitionException::what() const throw() {
    return _message.c_str();
}
