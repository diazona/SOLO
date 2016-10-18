#pragma once

#include <string>
#include <vector>
#include "../configuration/configuration.h"

/**
 * Stores high-level program configuration variables, e.g. information
 * about which command line options were passed.
 */
class ProgramConfiguration {
public:
                                      // char const * const * is one way to write the proper incantation to show that
                                      // this method won't modify the command-line arguments it gets passed
    ProgramConfiguration(const int argc, char const * const * argv);
    ~ProgramConfiguration();

    /** Indicates whether the --print-config or --no-print-config option was specified */
    bool print_config() const { return m_print_config; }
    /** Indicates whether the --print-integration-progress or --no-print-integration-progress option was specified */
    bool print_integration_progress() const { return m_print_integration_progress; }
    /** Indicates whether the --trace option was specified */
    bool trace() const { return m_trace; }
    /** Indicates whether the --trace-gdist option was specified */
    bool trace_gdist() const { return m_trace_gdist; }
    /** Indicates whether the --minmax option was specified */
    bool minmax() const { return m_minmax; }
    /** Indicates whether the --separate option was specified */
    bool separate() const { return m_separate; }

    double xg_min() const { return m_xg_min; }
    double xg_max() const { return m_xg_max; }

    const Configuration& config() const { return m_conf; }
private:
    /** Indicates whether the --print-config or --no-print-config option was specified */
    bool m_print_config;
    /** Indicates whether the --print-integration-progress or --no-print-integration-progress option was specified */
    bool m_print_integration_progress;
    /** Indicates whether the --trace option was specified */
    bool m_trace;
    /** Indicates whether the --trace-gdist option was specified */
    bool m_trace_gdist;
    /** Indicates whether the --minmax option was specified */
    bool m_minmax;
    /** Indicates whether the --separate option was specified */
    bool m_separate;
    /**
     * The configuration parameters to be used in the calculation. Information
     * collected from the command line options and read from configuration files
     * specified on the command line will be stored in this.
     */
    Configuration m_conf;

    double m_xg_min, m_xg_max;
};
