#pragma once

#include <ostream>
#include <string>
#include <vector>
#include "../configuration/context.h"
#include "../hardfactors/hardfactor.h"
#include "programconfiguration.h"

/**
 * Stores the results of the integration and contains methods to run the calculation.
 */
class ResultsCalculator {
public:
    // cc needs to be before tlctx and result_array_len because of initializer dependencies
    /** Collection of the contexts to be used for the calculation */
    const ContextCollection cc;

private:
    /** The thread-local context to be used for the calculation */
    ThreadLocalContext tlctx;

    /**
     * Stores all the hard factors and groups
     */
    HardFactorRegistry registry;

    /**
     * The list of groups of hard factors
     */
    std::vector<const HardFactorGroup*> hfgroups;

    /**
     * The list of names of the hard factors. They are stored according to the
     * group they were given in, then by order within the group.
     */
    std::vector<string> hfnames;

    /** The number of hard factor groups */
    size_t _hfglen;
    /** The number of hard factors */
    size_t _hflen;
    /** The length of the results arrays */
    size_t result_array_len;
    /** Flags the indices of results which have been successfully computed so far */
    bool* _valid;
    /** Array to hold the real parts of the results */
    double* real;
    /** Array to hold the imaginary parts of the results */
    double* imag;
    /** Array to hold the error bounds of the results */
    double* error;

    friend std::ostream& operator<<(std::ostream&, ResultsCalculator&);
public:
    /** Whether to trace execution */
    const bool trace;
    /** Whether to store minimum and maximum values */
    const bool minmax;
    /** Whether to calculate individual hard factors separately */
    const bool separate;
    /** Whether to print integration progress updates */
    const bool print_integration_progress;

    ResultsCalculator(const ProgramConfiguration& pc);
    ~ResultsCalculator();

    /**
     * Turns a context index and a hard factor group index into an index into a
     * 1D row-major array
     */
    size_t index_from(size_t ccindex, size_t hfindex);
    /**
     * Return whether the given combination of context index and hard factor
     * group index is valid - that is, whether a result has been computed
     * for that combination
     */
    bool valid(size_t ccindex, size_t hfindex);
    /**
     * Places the result at the given context index and hard factor group
     * index into the variables real, imag, and error. This should only
     * be called after calculate().
     */
    void result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error);
    /**
     * Runs the calculation.
     */
    void calculate();
private:
    /**
     * Parse the hard factor specifications collected in the constructor.
     *
     * This will force a ContextCollection to be constructed if it has not been
     * done already.
     */
    void parse_hf_specs(const std::vector<std::string>& hfspecs);

    /**
     * Construct an Integrator and use it
     */
    void integrate_hard_factor(const Context& ctx, const ThreadLocalContext& tlctx, const HardFactorList& hflist, size_t index);

    double xg_min, xg_max;
};
