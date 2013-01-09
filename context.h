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

#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "coupling.h"
#include "gluondist.h"

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
    double A;
    /** The centrality coefficient */
    double c;
    /** The exponent in the saturation scale formula */
    double lambda;
    /** The factorization scale */
    double mu2;
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
    
    /** GSL pseudorandom number generator algorithm */
    const gsl_rng_type* pseudorandom_generator_type;
    /** GSL pseudorandom number generator seed */
    unsigned long int pseudorandom_generator_seed;

    /** The gluon distribution */
    GluonDistribution* gdist;
    /** The coupling */
    Coupling* cpl;
    
    /** Number of MISER iterations
     * (unused unless integration strategy is MISER) */
    size_t miser_iterations;
    /** Number of VEGAS iterations when tuning the grid
     * (unused unless integration strategy is VEGAS) */
    size_t vegas_initial_iterations;
    /** Number of VEGAS iterations when actually integrating
     * (unused unless integration strategy is VEGAS) */
    size_t vegas_incremental_iterations;

    /** The precomputed value of c A^(1/3) Q_0^2 x_0^(lambda) */
    double Q02x0lambda;
    /** The precomputed value of tau = pT / sqs * exp(Y) */
    double tau;

    Context(
        double x0,
        double A,
        double c,
        double lambda,
        double mu2,
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
        size_t miser_iterations,
        size_t vegas_initial_iterations,
        size_t vegas_incremental_iterations,
        const gsl_rng_type* pseudorandom_generator_type,
        unsigned long int pseudorandom_generator_seed,
        GluonDistribution* gdist,
        Coupling* cpl) :
     x0(x0),
     A(A),
     c(c),
     lambda(lambda),
     mu2(mu2),
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
     miser_iterations(miser_iterations),
     vegas_initial_iterations(vegas_initial_iterations),
     vegas_incremental_iterations(vegas_incremental_iterations),
     pseudorandom_generator_type(pseudorandom_generator_type),
     pseudorandom_generator_seed(pseudorandom_generator_seed),
     gdist(gdist),
     cpl(cpl),
     Q02x0lambda(c * pow(A, 1.0d/3.0d) * pow(x0, lambda)),
     tau(sqrt(pT2)/sqs*exp(Y)) {}
    Context(const Context& other) :
     x0(other.x0),
     A(other.A),
     c(other.c),
     lambda(other.lambda),
     mu2(other.mu2),
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
     miser_iterations(other.miser_iterations),
     vegas_initial_iterations(other.vegas_initial_iterations),
     vegas_incremental_iterations(other.vegas_incremental_iterations),
     pseudorandom_generator_type(other.pseudorandom_generator_type),
     pseudorandom_generator_seed(other.pseudorandom_generator_seed),
     gdist(other.gdist),
     cpl(other.cpl),
     Q02x0lambda(other.Q02x0lambda),
     tau(other.tau) {}
};

/**
 * Write the Context out to the given output stream as a set of 
 *  key = value
 * lines. This is meant to be human-readable output.
 */
std::ostream operator<<(std::ostream out, Context& ctx);

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
      cpl(NULL) {
        setup_defaults();
    }
    /**
     * Construct a ContextCollection and initialize it with settings
     * read from the named file.
     */
    ContextCollection(const std::string& filename) :
      gdist(NULL),
      cpl(NULL) {
        setup_defaults();
        ifstream in(filename.c_str());
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
     * The gluon distribution. NULL until contexts are created.
     */
    GluonDistribution* gdist;
    /**
     * The coupling. NULL until contexts are created.
     */
    Coupling* cpl;

    friend std::istream& operator>>(std::istream& in, ContextCollection& cc);
    friend std::ostream& operator<<(std::ostream& out, ContextCollection& cc);
    friend class ThreadLocalContext;
    
private:
    /**
     * The map of key-value pairs provided to the ContextCollection.
     */
    std::multimap<std::string, std::string> options;
    /**
     * The Contexts. This is empty until contexts are created.
     */
    std::vector<Context> contexts;
    /**
     * Called from the constructor to set default values.
     */
    void setup_defaults();
    /**
     * Read a config file, or something in an equivalent format, from an input stream
     * and add the settings to the ContextCollection.
     */
    void read_config(std::istream& in);
    /**
     * Create the Context objects.
     */
    void create_contexts();
};

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
      pdf_object(new c_mstwpdf(cc.options.find("pdf_filename")->second.c_str())),
      ff_object(new DSSpiNLO(cc.options.find("ff_filename")->second.c_str())) {};
    
    ~ThreadLocalContext() {
        delete pdf_object;
        delete ff_object;
    }
};

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
        s << "No value for " << property << "!" << std::endl;
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
 * parsed into the correct type of object to pass to the Context constructor.
 */
template<typename T>
class InvalidPropertyValueException : public exception {
private:
    string _property;
    T _value;
    string _message;
public:
    InvalidPropertyValueException(const char* property, T& value) throw() : _property(property), _value(value) {
        std::ostringstream s;
        s << "Invalid value '" << value << "' for " << property << "!" << std::endl;
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

#endif // _CONTEXT_H_
