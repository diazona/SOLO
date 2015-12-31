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
#include <iostream>
#include <sstream>
#include <muParser.h>
#include <gsl/gsl_math.h>
#include "gsl_mu.h"
#include "hardfactor.h"
#include "hardfactor_parser.h"
#include "../typedefs.h"
#include "../utils/utils.h"

using mu::Parser;
using mu::value_type;
using mu::valmap_type;
using mu::varmap_type;
using std::cerr;

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

void ParsedHardFactorTerm::init_parser(Parser& parser, varmap_type& all_used_variables, const string& aux_variables, const string& real_expr, const string& imag_expr, const char* debug_message) {
    parser.SetExpr(aux_variables + e0(real_expr) + "," + e0(imag_expr));
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
    all_used_variables.insert(variables_from_parser.begin(), variables_from_parser.end());
}

ParsedHardFactorTerm::ParsedHardFactorTerm(
    const string& name,
    const string& implementation,
    const HardFactor::HardFactorOrder order,
    const IntegrationType* type,
    const string& Fs_real, const string& Fs_imag,
    const string& Fn_real, const string& Fn_imag,
    const string& Fd_real, const string& Fd_imag,
    const list<pair<string, string> >& variable_list) :
  m_name(name),
  m_implementation(implementation),
  m_order(order),
  mp_type(type) {
    aux_variable_count = variable_list.size();
    ostringstream aux_variables_oss;
    for (list<pair<string, string> >::const_iterator vi = variable_list.begin(); vi != variable_list.end(); vi++) {
        string var = vi->first;
        string expr = vi->second;
        aux_variables_oss << var << "=" << expr << ",";
    }
    string aux_variables = aux_variables_oss.str();
    varmap_type all_used_variables;
    init_parser(Fs_parser, all_used_variables, aux_variables, Fs_real, Fs_imag, "Fs");
    init_parser(Fn_parser, all_used_variables, aux_variables, Fn_real, Fn_imag, "Fn");
    init_parser(Fd_parser, all_used_variables, aux_variables, Fd_real, Fd_imag, "Fd");

    // make sure only existing variables are used
#define process(var) all_used_variables.erase(#var);
#include "../integration/ictx_var_list.inc"
#include "../configuration/ctx_var_list.inc"
#undef process
    aux_variable_storage = new double[aux_variable_count];
    size_t i = 0;
    for (list<pair<string, string> >::const_iterator vi = variable_list.begin(); vi != variable_list.end(); vi++, i++) {
        const string& varname = vi->first;
        Fs_parser.DefineVar(varname, &aux_variable_storage[i]);
        Fn_parser.DefineVar(varname, &aux_variable_storage[i]);
        Fd_parser.DefineVar(varname, &aux_variable_storage[i]);
        all_used_variables.erase(varname);
    }
    if (!all_used_variables.empty()) {
        ostringstream doss;
        doss << "Invalid variables used:";
        for (varmap_type::const_iterator it = all_used_variables.begin(); it != all_used_variables.end(); it++) {
            doss << " " << it->first;
        }
        delete[] aux_variable_storage;
        throw mu::ParserError(doss.str());
    }
}

ParsedHardFactorTerm::~ParsedHardFactorTerm() {
    delete[] aux_variable_storage;
}

// need to pull some tricks with templates to only call DefineVar on double-type variables
// see https://stackoverflow.com/a/16924234/56541
template<typename T> void define_one_variable(Parser& parser, const char* var_name, const T* var) {
}

template<> void define_one_variable(Parser& parser, const char* var_name, const double* var) {
    parser.DefineVar(var_name, const_cast<value_type*>(var));
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
#define process(var) define_one_variable(parser, #var, &(ictx->var));
#include "../integration/ictx_var_list.inc"
#undef process
#define process(var) define_one_variable(parser, #var, &(ictx->ctx->var));
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
    if (number_of_values < 2) {
        throw mu::ParserError("invalid number of values");
        return;
    }
    *real = values[number_of_values - 2];
    *imag = values[number_of_values - 1];
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

HardFactorParser::HardFactorParser(HardFactorRegistry& registry) :
    order(sentinel),
    type(NULL),
    registry(registry),
    hard_factor_callback(NULL),
    hard_factor_group_callback(NULL),
    error_handler(NULL) {
}

HardFactorParser::~HardFactorParser() {
    try {
        create_hard_factor_term();
    }
    catch (const IncompleteHardFactorDefinitionException& e) {
        // nothing we can do about it
        cerr << "Incomplete hard factor definition on destruction of parser" << endl;
    }
    reset_current_term();
    flush_groups();
}

string default_implementation(const IntegrationType* type) {
    string implementation;
    if (type == &position::dipole || type == &position::quadrupole) {
        implementation = "p";
    }
    else if (type == &radial::dipole || type == &radial::quadrupole || type == &radial::rescaled_dipole || type == &radial::rescaled_quadrupole) {
        implementation = "r";
    }
    else if (type == &momentum::momentum1 || type == &momentum::momentum2 || type == &momentum::momentum3 || type == &momentum::momentumxip1 ||
             type == &momentum::momentumxip2 || type == &momentum::none || type == &momentum::qlim || type == &momentum::radialmomentum1 ||
             type == &momentum::radialmomentum2 || type == &momentum::radialmomentum3) {
        implementation = "m";
    }
    assert(!implementation.empty());
    return implementation;
}

void HardFactorParser::parse_line(const string& line) {
    // an empty line signals the end of a hard factor term definition
    if (line.length() == 0) {
        create_hard_factor_term();
        return;
    }
    // ignore comments
    if (line[0] == '#') {
        return;
    }
    /* A hard factor definition consists of a name, an order
     * (LO or NLO), an integration type, and six functions
     * giving the real and imaginary components of Fs, Fn, Fd.
     * Each part of the definition is given as name = value.
     */
    vector<string> parts = split(line, "=", 2);
    assert(parts.size() > 0);
    if (parts.size() == 1) {
        unparsed_hard_factor_group_specs.push_back(parts[0]);
        return;
    }
    else if (parts.size() != 2) {
        throw InvalidHardFactorDefinitionException(line, parts[0], parts.size() > 1 ? parts[1] : "", "Incomplete or malformed property");
    }
    string key = trim(parts[0]);
    string value = trim(parts[1]);
    if (key == "name") {
        split_hardfactor(value, name, implementation);
    }
    else if (key == "order") {
        if (value == "lo") {
            order = HardFactor::LO;
        }
        else if (value == "nlo") {
            order = HardFactor::NLO;
        }
        /* no support for mixed-order hard factors because
         * this function only creates HardFactorTerms which
         * can only be LO or NLO */
        else if (value == "mixed") {
            throw InvalidHardFactorDefinitionException(line, key, value, "Mixed-order hard factors not supported:");
        }
        else {
            throw InvalidHardFactorDefinitionException(line, key, value);
        }
    }
    else if (key == "type") {
        if (value == "none") {
            type = &momentum::none;
        }
        else if (value == "momentum1") {
            type = &momentum::momentum1;
        }
        else if (value == "momentum2") {
            type = &momentum::momentum2;
        }
        else if (value == "momentum3") {
            type = &momentum::momentum3;
        }
        else if (value == "radial momentum1") {
            type = &momentum::radialmomentum1;
        }
        else if (value == "radial momentum2") {
            type = &momentum::radialmomentum2;
        }
        else if (value == "radial momentum3") {
            type = &momentum::radialmomentum3;
        }
        else if (value == "momentumxip1") {
            type = &momentum::momentumxip1;
        }
        else if (value == "momentumxip2") {
            type = &momentum::momentumxip2;
        }
        else if (value == "qlim") {
            type = &momentum::qlim;
        }
        else if (value == "cartesian dipole") {
            type = &position::dipole;
        }
        else if (value == "cartesian quadrupole") {
            type = &position::quadrupole;
        }
        else if (value == "radial dipole") {
            type = &radial::dipole;
        }
        else if (value == "radial quadrupole") {
            type = &radial::quadrupole;
        }
        else if (value == "radial rescaled dipole") {
            type = &radial::rescaled_dipole;
        }
        else if (value == "radial rescaled quadrupole") {
            type = &radial::rescaled_quadrupole;
        }
        else {
            throw InvalidHardFactorDefinitionException(line, key, value, "Unknown integration type:");
        }
    }
    else if (key == "Fs real") {
        Fs_real = value;
    }
    else if (key == "Fs imag") {
        Fs_imag = value;
    }
    else if (key == "Fn real") {
        Fn_real = value;
    }
    else if (key == "Fn imag") {
        Fn_imag = value;
    }
    else if (key == "Fd real") {
        Fd_real = value;
    }
    else if (key == "Fd imag") {
        Fd_imag = value;
    }
    /* Anything other than those keys, if it doesn't contain spaces, is
     * taken to represent either:
     *
     * - A definition of a variable to be used in the hard factor expressions,
     *   if it occurs when a hard factor term definition is already in progress
     *   (i.e. any of name, order, type, etc. have already been given but the
     *   hard factor definition is not complete); or
     * - A specification of a hard factor which contains multiple terms,
     *   if there is not a hard factor definition already in progress.
     */
    else if (key.find_first_of(" \t\r\n") == string::npos) {
        if (hard_factor_definition_empty()) {
            parse_composite_hard_factor(key, value);
        }
        else {
            parse_variable_definition(key, value);
        }
    }
    else {
        throw InvalidHardFactorDefinitionException(line, key, value, "Unknown property:");
    }
}

void HardFactorParser::parse_file(const string& filename) {
    ifstream hfinput(filename.c_str());
    string line;
    size_t i;

    for (i = 0, getline(hfinput, line); hfinput; i++, getline(hfinput, line)) {
        try {
            parse_line(line);
        }
        catch (const InvalidHardFactorDefinitionException& e) {
            if (error_handler == NULL || error_handler(e, filename, i)) {
                throw;
            }
        }
        catch (const IncompleteHardFactorDefinitionException& e) {
            if (error_handler == NULL || error_handler(e, filename, i)) {
                throw;
            }
        }
    }
    create_hard_factor_term();
}

void HardFactorParser::parse_variable_definition(const string& key, const string& value) {
    pair<string, string> p(key, value);
    variable_definitions.push_back(p);
}


const ParsedCompositeHardFactor* HardFactorParser::parse_composite_hard_factor(const string& key, const string& value) {
    assert(hard_factor_definition_empty());
    vector<string> term_labels = split(value, ",+");
    if (term_labels.empty()) {
        throw InvalidHardFactorSpecException(key, "no hard factor terms provided in definition");
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
        string name, implementation;
        split_hardfactor(s, name, implementation);
        const HardFactor* hf;
        if (implementation.empty()) {
            hf = registry.get_hard_factor(name);
        }
        else {
            hf = registry.get_hard_factor(name, implementation);
        }
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
        throw InvalidHardFactorSpecException(value, "no valid hard factors provided in definition");
    }

    split_hardfactor(key, name, implementation);
    ParsedCompositeHardFactor* hf = new ParsedCompositeHardFactor(name, order, hftlist);

    if (implementation.empty()) {
        // temporary hack: take the type of the composite hard factor
        // to be the type of its first term
        implementation = default_implementation(hftlist[0]->get_type());
    }

    registry.add_hard_factor(name, implementation, hf, true);
    hard_factors.push_back(hf);
    if (hard_factor_callback) {
        hard_factor_callback(*hf);
    }
    reset_current_term();
    return hf;
}


const HardFactorGroup* HardFactorParser::parse_hard_factor_group(const string& spec) {
    string specname(spec);
    string specbody(spec);

    vector<string> namesplitspec = split(spec, ":", 2);
    assert(namesplitspec.size() == 1 || namesplitspec.size() == 2);
    if (namesplitspec.size() > 1) {
        if (namesplitspec[0].size() == 0) {
            // this means there was simply a leading colon; ignore it
            specname = namesplitspec[1];
        }
        else {
            specname = namesplitspec[0];
        }
        specbody = namesplitspec[1];
    }

    vector<string> splitspec;
    splitspec = split(specbody, ", ");

    HardFactorList hfobjs;
    // Iterate over the individual hard factor names
    for (vector<string>::iterator it = splitspec.begin(); it != splitspec.end(); it++) {
        string name, implementation;
        split_hardfactor(*it, name, implementation);
        const HardFactor* hf;
        if (implementation.empty()) {
            hf = registry.get_hard_factor(name);
        }
        else {
            hf = registry.get_hard_factor(name, implementation);
        }
        if (hf == NULL) {
            // the string failed to parse
            throw InvalidHardFactorSpecException(spec, "No such hard factor");
        }
        else {
            hfobjs.push_back(hf);
        }
    }
    if (hfobjs.size() == 0) {
        throw InvalidHardFactorSpecException(spec, "No valid hard factors in specification");
    }
    else {
        // parsing seems to have succeeded, so go ahead and create the object
        return new HardFactorGroup(specname, hfobjs, splitspec);
    }
}

void HardFactorParser::flush_groups() {
    while (!unparsed_hard_factor_group_specs.empty()) {
        string line = unparsed_hard_factor_group_specs.front();
        const HardFactorGroup* hfg = parse_hard_factor_group(line);
        if (hfg != NULL) {
            if (registry.get_hard_factor_group(hfg->label) != NULL) {
                throw InvalidHardFactorDefinitionException(line, hfg->label, line, "Duplicate hard factor group label");
            }
            registry.add_hard_factor_group(hfg, true);
            if (hard_factor_group_callback) {
                hard_factor_group_callback(*hfg);
            }
        }
        else {
            throw InvalidHardFactorSpecException(line, "Unknown error parsing group specification (parser returned NULL)");
        }
        unparsed_hard_factor_group_specs.pop_front();
    }
}


void HardFactorParser::split_hardfactor(const string& spec, string& name, string& implementation) {
    vector<string> splitspec;
    splitspec = split(spec, ".", 2);
    assert(splitspec.size() == 1 || splitspec.size() == 2);
    name = splitspec[0];
    if (splitspec.size() == 2) {
        implementation = splitspec[1];
    }
    else {
        implementation = "";
    }
}

void HardFactorParser::set_error_handler(bool (*error_handler)(const exception&, const string&, const size_t)) {
    this->error_handler = error_handler;
}

void HardFactorParser::set_hard_factor_callback(void (*callback)(const HardFactor&)) {
    this->hard_factor_callback = callback;
}

void HardFactorParser::set_hard_factor_group_callback(void (*callback)(const HardFactorGroup&)) {
    this->hard_factor_group_callback = callback;
}


bool HardFactorParser::hard_factor_definition_empty() const {
    return
      type == NULL &&
      order == sentinel &&
      name.empty() &&
      implementation.empty() &&
      Fs_real.empty() && Fs_imag.empty() &&
      Fn_real.empty() && Fn_imag.empty() &&
      Fd_real.empty() && Fd_imag.empty() &&
      variable_definitions.empty();
}

bool HardFactorParser::hard_factor_definition_complete() const {
    return
      type != NULL &&
      order != sentinel &&
      !name.empty();
}

const ParsedHardFactorTerm* HardFactorParser::create_hard_factor_term() {
    if (hard_factor_definition_empty()) {
        return NULL;
    }
    if (!hard_factor_definition_complete()) {
        throw IncompleteHardFactorDefinitionException();
    }

    ParsedHardFactorTerm* hf = new ParsedHardFactorTerm(name, implementation.empty() ? default_implementation(type) : implementation, order, type, Fs_real, Fs_imag, Fn_real, Fn_imag, Fd_real, Fd_imag, variable_definitions);
    registry.add_hard_factor(hf, true);
    hard_factors.push_back(hf);
    if (hard_factor_callback != NULL) {
        hard_factor_callback(*hf);
    }
    reset_current_term();
    return hf;
}

void HardFactorParser::reset_current_term() {
    name.erase();
    implementation.erase();
    order = sentinel;
    type = NULL;
    variable_definitions.clear();
    Fs_real.clear();
    Fs_imag.clear();
    Fn_real.clear();
    Fn_imag.clear();
    Fd_real.clear();
    Fd_imag.clear();
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

InvalidHardFactorSpecException::InvalidHardFactorSpecException(const string& hfspec, const string& message) throw() : hfspec(hfspec) {
    std::ostringstream s;
    s << message << "; spec " << hfspec;
    _message = s.str();
}

void InvalidHardFactorSpecException::operator=(const InvalidHardFactorSpecException& other) {
    _message = other._message;
}

const char* InvalidHardFactorSpecException::what() const throw() {
    return _message.c_str();
}
