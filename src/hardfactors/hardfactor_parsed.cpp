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
#include "gsl_mu.h"
#include "hardfactor.h"
#include "hardfactors_position.h"
#include "hardfactors_radial.h"
#include "hardfactors_momentum.h"
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

ParsedHardFactorTerm::ParsedHardFactorTerm(
    const string& name, const HardFactor::HardFactorOrder order, const IntegrationType* type, const string& Fs_real, const string& Fs_imag, const string& Fn_real, const string& Fn_imag, const string& Fd_real, const string& Fd_imag) :
  m_name(name),
  m_order(order),
  mp_type(type) {
    Fs_parser.SetExpr(e0(Fs_real) + "," + e0(Fs_imag));
    Fn_parser.SetExpr(e0(Fn_real) + "," + e0(Fn_imag));
    Fd_parser.SetExpr(e0(Fd_real) + "," + e0(Fd_imag));
    
    mu_load_gsl(Fs_parser);
    mu_load_gsl(Fn_parser);
    mu_load_gsl(Fd_parser);
    
    // use GetUsedVar() to trigger an attempt to parse the expression
    // we want any parsing errors to be revealed now, not when evaluating expressions
    varmap_type all_variables;
    varmap_type variables_from_parser;
    variables_from_parser = Fs_parser.GetUsedVar();
    all_variables.insert(variables_from_parser.begin(), variables_from_parser.end());
    variables_from_parser = Fn_parser.GetUsedVar();
    all_variables.insert(variables_from_parser.begin(), variables_from_parser.end());
    variables_from_parser = Fd_parser.GetUsedVar();
    all_variables.insert(variables_from_parser.begin(), variables_from_parser.end());
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
}

void ParsedHardFactorTerm::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx != m_ictx_Fn) {
        m_ictx_Fn = ictx;
        define_variables(Fn_parser, ictx);
    }
    evaluate_hard_factor(Fn_parser, real, imag);
}

void ParsedHardFactorTerm::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx != m_ictx_Fd) {
        m_ictx_Fd = ictx;
        define_variables(Fd_parser, ictx);
    }
    evaluate_hard_factor(Fd_parser, real, imag);
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
        vector<string> term_labels = split(parts[1], ",");
        if (term_labels.empty()) {
            throw InvalidHardFactorSpecException(parts[0], "no hard factor terms provided in definition");
        }
        HardFactorTermList hftlist;
        HardFactor::HardFactorOrder order = sentinel;
        for (vector<string>::const_iterator it = term_labels.begin(); it != term_labels.end(); it++) {
            if (it->empty()) {
                continue;
            }
            // possible future enhancement: handle the case where *it names
            // a composite hard factor by extracting the terms and adding them
            // individually
            const HardFactor* hf = parse_hardfactor(*it);
            if (hf == NULL) {
                throw InvalidHardFactorSpecException(*it, "no hard factor with that name");
            }
            const HardFactorTerm* hft = dynamic_cast<const HardFactorTerm*>(hf);
            if (hft == NULL) {
                throw InvalidHardFactorSpecException(*it, "can't handle composite hard factors");
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
    
    if (hard_factor_definition_complete()) {
        create_hard_factor_term();
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
                throw e;
            }
        }
        catch (const IncompleteHardFactorDefinitionException& e) {
            if (error_handler == NULL || error_handler(e, line, i)) {
                throw e;
            }
        }
    }
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
}

