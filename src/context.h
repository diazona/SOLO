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

#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include <cassert>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "coupling.h"
#include "factorizationscale.h"
#include "gluondist.h"

/**
 * Enumerates the types of Monte Carlo integration available
 */
typedef enum {MC_PLAIN, MC_MISER, MC_VEGAS, MC_QUASI} integration_strategy;

/**
 * Enumerates the types of projectile available
 */
typedef enum {proton, deuteron} projectile_type;

/**
 * An exception to throw when creating a Context, if the Context constructor
 * requires a property that hasn't been added to the ContextCollection and
 * has no default value.
 */
class MissingPropertyException : public exception {
private:
    string _property;
    string _message;
public:
    MissingPropertyException(const char* property) throw() : _property(property) {
        std::ostringstream s;
        s << "No value for " << property << "!";
        _message = s.str();
    }
    MissingPropertyException(const MissingPropertyException& other) throw() : _property(other._property), _message(other._message) {}
    ~MissingPropertyException() throw() {}
    void operator=(const MissingPropertyException& other) {
        _property = other._property;
        _message = other._message;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * An exception to throw when creating a Context, if the value can't be
 * parsed into the correct type of object to pass to the Context constructor
 * or if the parsed value is invalid for that property for some other reason.
 *
 * If it's not clear from the string representation of the value why it's
 * inappropriate for the propery, an extra message should be added clarifying
 * that.
 */
template<typename T>
class InvalidPropertyValueException : public exception {
private:
    const string _property;
    const T _value;
    string _message;
public:
    InvalidPropertyValueException(const char* property, const T& value) throw() : _property(property), _value(value) {
        std::ostringstream s;
        s << "Invalid value '" << value << "' for " << property << "!";
        _message = s.str();
    }
    InvalidPropertyValueException(const char* property, const T& value, const char* extra_message) throw() : _property(property), _value(value) {
        std::ostringstream s;
        s << "Invalid value '" << value << "' for " << property << "!" << extra_message;
        _message = s.str();
    }
    InvalidPropertyValueException(const InvalidPropertyValueException& other) throw() : _property(other._property), _value(other._value), _message(other._message) {}
    ~InvalidPropertyValueException() throw() {}
    void operator=(const InvalidPropertyValueException& other) throw() {
        _property = other._property;
        _value = other._value;
        _message = other._message;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * An exception to throw when creating a Context, if the values of the kinematic
 * variables are physically invalid or inconsistent.
 */
class InvalidKinematicsException : public exception {
private:
    string _message;
public:
    InvalidKinematicsException(const char* message) throw() : _message(message) {
    }
    InvalidKinematicsException(const InvalidKinematicsException& other) throw() : _message(other._message) {}
    ~InvalidKinematicsException() throw() {}
    void operator=(const InvalidKinematicsException& other) {
        _message = other._message;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * Storage for all the assorted parameters that get used in the integration.
 * 
 * There is a different Context object for each distinct set of parameters
 * that a result is computed for.
 * 
 * Other code should treat the contents of a Context as constant, and not
 * modify it.
 */
class Context {
public:
    // Actually making these const is tricky
    /** The fit parameter from the saturation scale */
    double x0;
    /** The mass number */
    double mass_number;
    /** The centrality coefficient */
    double centrality;
    /** The exponent in the saturation scale formula */
    double lambda;
    /** The number of colors */
    double Nc;
    /** The number of flavors */
    double Nf;
    /** The color factor */
    double CF;
    double TR;
    /** The nucleon cross sectional area */
    double Sperp;
    /** The transverse momentum squared */
    double pT2;
    /** The collider CM energy */
    double sqs;
    /** The rapidity */
    double Y;
    /** Name of the file PDF data was read from */
    std::string pdf_filename;
    /** Name of the file FF data was read from */
    std::string ff_filename;
    
    /** GSL quasirandom number generator algorithm */
    const gsl_qrng_type* quasirandom_generator_type;
    /** GSL pseudorandom number generator algorithm */
    const gsl_rng_type* pseudorandom_generator_type;
    /** GSL pseudorandom number generator seed */
    unsigned long int pseudorandom_generator_seed;

    /** The gluon distribution */
    GluonDistribution* gdist;
    /** The coupling */
    Coupling* cpl;
    /** The factorization scale */
    FactorizationScale* fs;
    /** Whether to apply the optimization that sets ln(c_0^2/(r^2 mu^2)) to zero */
    bool c0r_optimization;
    /**
     * Whether to use the prescription from Collins/Soper/Sterman for modifying r,
     * from DOI 10.1016/0550-3213(85)90479-1
     */
    bool css_r_regularization;
    /** The cutoff value for the CSS r regularization */
    double css_r2_max;
    
    /** The factor in front of the resummation term, see H1qqCorrection */
    double resummation_constant;
    
    /** Whether to use exact (or approximate) kinematic expressions */
    bool exact_kinematics;

    /** Projectile type */
    projectile_type projectile;
    /** Product hadron */
    DSSpiNLO::hadron hadron;
    /** The type of integration to be used */
    integration_strategy strategy;
    
    /** Maximum allowed absolute error, for integration strategies that use it */
    double abserr;
    /** Maximum allowed relative error, for integration strategies that use it */
    double relerr;
    
    /** Number of iterations for cubature */
    size_t cubature_iterations;
    /** Number of MISER iterations
     * (unused unless integration strategy is MISER) */
    size_t miser_iterations;
    /** Number of VEGAS iterations when tuning the grid
     * (unused unless integration strategy is VEGAS) */
    size_t vegas_initial_iterations;
    /** Number of VEGAS iterations when actually integrating
     * (unused unless integration strategy is VEGAS) */
    size_t vegas_incremental_iterations;
    /** Number of iterations in quasi Monte Carlo
     * (unusued unless integration strategy is QUASI) */
    size_t quasi_iterations;

    /** The limit of integration over infinite regions */
    double inf;

    /** The precomputed value of c A^(1/3) Q_0^2 x_0^(lambda) */
    double Q02x0lambda;
    /** The precomputed value of tau = pT / sqs * exp(Y) */
    double tau;
    /** The precomputed value of tauhat = pT / sqs * (exp(Y) + exp(-Y)) */
    double tauhat;
    
    Context(
        double x0,
        double mass_number,
        double centrality,
        double lambda,
        double Nc,
        double Nf,
        double CF,
        double TR,
        double Sperp,
        double pT2,
        double sqs,
        double Y,
        std::string pdf_filename,
        std::string ff_filename,
        projectile_type projectile,
        DSSpiNLO::hadron hadron,
        integration_strategy strategy,
        size_t cubature_iterations,
        size_t miser_iterations,
        size_t vegas_initial_iterations,
        size_t vegas_incremental_iterations,
        size_t quasi_iterations,
        double abserr,
        double relerr,
        const gsl_qrng_type* quasirandom_generator_type,
        const gsl_rng_type* pseudorandom_generator_type,
        unsigned long int pseudorandom_generator_seed,
        GluonDistribution* gdist,
        Coupling* cpl,
        FactorizationScale* fs,
        bool c0r_optimization,
        bool css_r_regularization,
        double css_r2_max,
        double resummation_constant,
        bool exact_kinematics,
        double inf) :
     x0(x0),
     mass_number(mass_number),
     centrality(centrality),
     lambda(lambda),
     Nc(Nc),
     Nf(Nf),
     CF(CF),
     TR(TR),
     Sperp(Sperp),
     pT2(pT2),
     sqs(sqs),
     Y(Y),
     pdf_filename(pdf_filename),
     ff_filename(ff_filename),
     projectile(projectile),
     hadron(hadron),
     strategy(strategy),
     cubature_iterations(cubature_iterations),
     miser_iterations(miser_iterations),
     vegas_initial_iterations(vegas_initial_iterations),
     vegas_incremental_iterations(vegas_incremental_iterations),
     quasi_iterations(quasi_iterations),
     abserr(abserr),
     relerr(relerr),
     quasirandom_generator_type(quasirandom_generator_type),
     pseudorandom_generator_type(pseudorandom_generator_type),
     pseudorandom_generator_seed(pseudorandom_generator_seed),
     gdist(gdist),
     cpl(cpl),
     fs(fs),
     c0r_optimization(c0r_optimization),
     css_r_regularization(css_r_regularization),
     css_r2_max(css_r2_max),
     resummation_constant(resummation_constant),
     exact_kinematics(exact_kinematics),
     inf(inf),
     Q02x0lambda(centrality * pow(mass_number, 1.0d/3.0d) * pow(x0, lambda)),
     tau(sqrt(pT2)/sqs*exp(Y)),
     tauhat(sqrt(pT2)/sqs*(exp(Y)+exp(-Y))) {
        if (tau > 1) {
            throw InvalidKinematicsException("τ > 1: empty phase space");
        }
        if (tauhat > 1) {
            throw InvalidKinematicsException("\\hat{τ} > 1: empty phase space");
        }
        if (css_r_regularization) {
            assert(css_r2_max > 0);
        }
    }
    Context(const Context& other) :
     x0(other.x0),
     mass_number(other.mass_number),
     centrality(other.centrality),
     lambda(other.lambda),
     Nc(other.Nc),
     Nf(other.Nf),
     CF(other.CF),
     TR(other.TR),
     Sperp(other.Sperp),
     pT2(other.pT2),
     sqs(other.sqs),
     Y(other.Y),
     pdf_filename(other.pdf_filename),
     ff_filename(other.ff_filename),
     projectile(other.projectile),
     hadron(other.hadron),
     strategy(other.strategy),
     cubature_iterations(other.cubature_iterations),
     miser_iterations(other.miser_iterations),
     vegas_initial_iterations(other.vegas_initial_iterations),
     vegas_incremental_iterations(other.vegas_incremental_iterations),
     quasi_iterations(other.quasi_iterations),
     abserr(other.abserr),
     relerr(other.relerr),
     quasirandom_generator_type(other.quasirandom_generator_type),
     pseudorandom_generator_type(other.pseudorandom_generator_type),
     pseudorandom_generator_seed(other.pseudorandom_generator_seed),
     gdist(other.gdist),
     cpl(other.cpl),
     fs(other.fs),
     c0r_optimization(other.c0r_optimization),
     css_r_regularization(other.css_r_regularization),
     css_r2_max(other.css_r2_max),
     resummation_constant(other.resummation_constant),
     exact_kinematics(other.exact_kinematics),
     inf(other.inf),
     Q02x0lambda(other.Q02x0lambda),
     tau(other.tau),
     tauhat(other.tauhat) {}
};

/**
 * The "context factory" and a repository for all settings.
 * 
 * A ContextCollection is able to read a configuration file in
 *  key = value
 * format, and store all the settings read. It allows multiple values
 * of pT and/or Y, but only one value of any other setting.
 * 
 * After all configuration files have been read, the ContextCollection
 * can be used to create a list of Context objects, one for each
 * combination of pT and Y. Calling any of the accessor methods
 * (get_context(), operator[](), begin(), end()) causes the set of
 * Contexts to be created, and also freezes the ContextCollection so
 * that the settings it holds can no longer be modified.
 * 
 * Several methods are named similar to, and behave similar to, their
 * counterparts in std::vector, allowing a ContextCollection to be
 * indexed or iterated over much like a vector. (It should be considered
 * read-only; do not assign to the iterator.)
 */
class ContextCollection {
public:
    typedef std::vector<Context>::iterator iterator;

    ContextCollection() :
      gdist(NULL),
      cpl(NULL),
      fs(NULL),
      trace_gdist(false),
      contexts_created(false) {
        setup_defaults();
    }
    /**
     * Construct a ContextCollection and initialize it with settings
     * read from the named file.
     */
    ContextCollection(const std::string& filename) :
      gdist(NULL),
      cpl(NULL),
      fs(NULL),
      trace_gdist(false),
      contexts_created(false) {
        setup_defaults();
        ifstream in(filename.c_str());
        if (!in.good()) {
            throw ios_base::failure("Unable to read file " + filename);
        }
        read_config(in);
        in.close();
    }
    ~ContextCollection() {
        delete gdist;
        gdist = NULL;
        delete cpl;
        cpl = NULL;
    }
    /**
     * Gets the nth context. Internally, the contexts are stored indexed
     * first by pT and then by Y. That is, if N is the number of Y values:
     * - Context 0 has pT[0] and Y[0]
     * - Context 1 has pT[0] and Y[1]
     * ...
     * - Context N-1 has pT[0] and Y[N-1]
     * - Context N has pT[1] and Y[0]
     * and so on.
     * 
     * When this method is called, if the Context objects have not already
     * been created, this creates the Contexts and freezes the ContextCollection.
     */
    Context& get_context(size_t n);
    /**
     * Allows access to contexts by subscript notation.
     * 
     * This is exactly equivalent to get_context().
     * 
     * When this method is called, if the Context objects have not already
     * been created, this creates the Contexts and freezes the ContextCollection.
     */
    Context& operator[](size_t n);
    /**
     * Tests whether the ContextCollection is empty.
     * 
     * This will return true if the number of pT values or the number
     * of Y values held by the ContextCollection is zero.
     * 
     * If this returns true, calling get_context() with any argument
     * will cause an error.
     */
    bool empty();
    /**
     * Returns the size of this ContextCollection.
     * 
     * This returns the product of the number of pT values specified so far
     * and the number of Y values specified so far. Before Contexts are created,
     * the return value can change as more settings are added. After Contexts
     * are created, the return value will not change.
     */
    size_t size();
    /**
     * Returns an iterator to the first Context.
     * 
     * When this method is called, if the Context objects have not already
     * been created, this creates the Contexts and freezes the ContextCollection.
     */
    iterator begin();
    /**
     * Returns an iterator to one past the last Context.
     * 
     * When this method is called, if the Context objects have not already
     * been created, this creates the Contexts and freezes the ContextCollection.
     */
    iterator end();
    
    /**
     * Removes all settings with the given key.
     */
    void erase(std::string key);
    /**
     * Returns the value of a setting, if any, or the empty string if unset.
     */
    std::string get(std::string key, size_t index = 0);
    /**
     * Add a setting, replacing any existing settings with the same key.
     */
    void set(std::string key, std::string value);
    /**
     * Add a setting. If the key is "pT" or "Y" (case insensitive), any existing
     * settings with the same key are left alone. For other keys, this behaves
     * identically to set().
     */
    void add(std::string key, std::string value);
    
    /**
     * Create the Context objects.
     */
    void create_contexts();

    /**
     * Read a config file, or something in an equivalent format, from an input stream
     * and add the settings to the ContextCollection.
     */
    void read_config(std::istream& in);
    /**
     * Process a string representing one line of a config file (i.e. one setting)
     */
    void read_config_line(std::string& line);
    
    /**
     * Whether to use the tracing gluon distribution wrapper. (See gluondist.h/cpp)
     * Changes made to this variable after contexts are created have no effect.
     */
    bool trace_gdist;

    friend ContextCollection& operator>>(std::string& line, ContextCollection& cc);
    friend std::istream& operator>>(std::istream& in, ContextCollection& cc);
    friend std::ostream& operator<<(std::ostream& out, ContextCollection& cc);
    friend class ThreadLocalContext;
    
private:
    /**
     * The gluon distribution. NULL until contexts are created.
     */
    GluonDistribution* gdist;
    /**
     * The coupling. NULL until contexts are created.
     */
    Coupling* cpl;
    /**
     * The factorization scale strategy. NULL until contexts are created.
     */
    FactorizationScale* fs;
    /**
     * The map of key-value pairs provided to the ContextCollection.
     */
    std::multimap<std::string, std::string> options;
    /**
     * The Contexts. This is empty until contexts are created.
     */
    std::vector<Context> contexts;
    /**
     * Whether contexts have been created.
     */
    bool contexts_created;
    /**
     * Called from the constructor to set default values.
     */
    void setup_defaults();
    
    /* Auxiliary methods and variables used to create gluon distributions */
    double Q02, x0, lambda, sqs, inf;
    vector<double> pT, Y;
    GBWGluonDistribution* create_gbw_gluon_distribution();
    MVGluonDistribution* create_mv_gluon_distribution();
    FixedSaturationMVGluonDistribution* create_fmv_gluon_distribution();
    PlateauPowerGluonDistribution* create_pp_gluon_distribution();
    FileDataGluonDistribution* create_file_gluon_distribution(GluonDistribution* lower_dist, GluonDistribution* upper_dist, const bool extended);
    GluonDistribution* create_gluon_distribution(const string&);
};

/**
 * Allows adding a setting to a ContextCollection using >> notation.
 */
ContextCollection& operator>>(std::string& line, ContextCollection& cc);
/**
 * Allows reading a ContextCollection in from a stream using >> notation.
 */
std::istream& operator>>(std::istream& in, ContextCollection& cc);
/**
 * Allows writing a ContextCollection out to a stream using << notation.
 * 
 * What is written out is just the list of key-value pairs. The output
 * could be read in to reconstruct the ContextCollection.
 */
std::ostream& operator<<(std::ostream& out, ContextCollection& cc);

/**
 * Allows writing a Context out to a stream using << notation.
 * 
 * This is a human-readable representation and cannot necessarily
 * be used to reconstruct the Context programmatically.
 */
std::ostream& operator<<(std::ostream& out, Context& ctx);

/**
 * Another Context-like class that holds objects which should not be shared
 * among threads or processes.
 * 
 * In the current state of the program, there isn't any particular reason
 * to have this, because there is no multithreading or multiprocessing being
 * used.
 */
class ThreadLocalContext {
private:
    friend class IntegrationContext;
    /**
     * The object that holds the PDF data
     */
    c_mstwpdf* pdf_object;
    /**
     * The object that holds the FF data
     */
    DSSpiNLO* ff_object;

public:
    ThreadLocalContext(const char* pdf_filename, const char* ff_filename) :
      pdf_object(new c_mstwpdf(pdf_filename)),
      ff_object(new DSSpiNLO(ff_filename)) {};
    ThreadLocalContext(const Context& ctx) :
      pdf_object(new c_mstwpdf(ctx.pdf_filename.c_str())),
      ff_object(new DSSpiNLO(ctx.ff_filename.c_str())) {};
    ThreadLocalContext(const ContextCollection& cc) :
      pdf_object(NULL), ff_object(NULL) {
        multimap<string,string>::const_iterator it = cc.options.find("pdf_filename");
        if (it == cc.options.end()) {
            throw "no PDF filename";
        }
        pdf_object = new c_mstwpdf(it->second.c_str());
        it = cc.options.find("ff_filename");
        if (it == cc.options.end()) {
            throw "no FF filename";
        }
        ff_object = new DSSpiNLO(it->second.c_str());
    };
    
    ~ThreadLocalContext() {
        delete pdf_object;
        delete ff_object;
    }
};


#endif // _CONTEXT_H_
