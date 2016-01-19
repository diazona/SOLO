#include <string>
#include <vector>
#include "../log.h"
#include "../configuration/context.h"
#include "../utils/utils.h"
#include "programconfiguration.h"
#include "trace.h"

using std::endl;
using std::cerr;
using std::string;
using std::vector;
using trace_variable::trace_vars;

ProgramConfiguration::ProgramConfiguration(const int argc, char const * const * argv) :
    m_trace(false),
    m_trace_gdist(false),
    m_minmax(false),
    m_separate(false),
    m_print_config(true),
    m_print_integration_progress(true),
    m_conf(canonicalize),
    m_xg_min(0),
    m_xg_max(1)
{
    string gdist_type;

    // parse command line options
    bool current_arg_is_config_line = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (current_arg_is_config_line) {
            m_conf.read_config_line(a);
            current_arg_is_config_line = false;
        }
        else if (a.compare(0, 10, "--ygrange=") == 0) {
            vector<string> v = split(a, "=", 2);
            if (v.size() == 2) {
                vector<string> r = split(v[1], ":", 2);
                if (r.size() == 2) {
                    m_xg_min = exp(-atof(r[1].c_str()));
                    m_xg_max = exp(-atof(r[0].c_str()));
                    if (m_xg_min > m_xg_max) {
                        cerr << "WARNING: reversing inverted range for ln(1/xg)" << endl;
                        double t = m_xg_min;
                        m_xg_min = m_xg_max;
                        m_xg_max = t;
                    }
                }
                else {
                    cerr << "invalid range for ln(1/xg): " << v[1] << endl;
                }
            }
        }
        else if (a.compare(0, 8, "--trace=") == 0) {
            vector<string> v = split(a, "=", 2);
            if (v.size() == 2) {
                if (v[1] == "*" || v[1] == "all") {
                    trace_vars.set();
                }
                else {
                    vector<string> w = split(v[1], ",");
                    for (vector<string>::iterator it = w.begin(); it != w.end(); it++) {
                        bool handled = false;
                        #define process(v) if (*it == #v) { trace_vars.set(static_cast<size_t>(trace_variable::v)); handled = true; }
                        #include "../integration/ictx_var_list.inc"
                        #undef process
                        if (!handled) {
                            cerr << "unknown trace variable " << *it << endl;
                        }
                    }
                }
                m_trace = trace_vars.any();
            }
        }
        else if (a == "--trace") {
            m_trace = true;
        }
        else if (a == "--trace-gdist") {
            m_trace_gdist = true;
        }
        else if (a == "--print-config") {
            m_print_config = true;
        }
        else if (a == "--no-print-config") {
            m_print_config = false;
        }
        else if (a == "--print-integration-progress") {
            m_print_integration_progress = true;
        }
        else if (a == "--no-print-integration-progress") {
            m_print_integration_progress = false;
        }
        else if (a == "--minmax") {
            m_minmax = true;
        }
        else if (a == "--separate") {
            m_separate = true;
        }
        else if (a == "-o" || a == "--option") {
            current_arg_is_config_line = true;
        }
        else if (a.compare(0, 2, "-o") == 0) {
            string s(a.substr(2));
            m_conf.read_config_line(s);
        }
        else if (a.compare(0, 8, "--option") == 0) {
            string s(a.substr(8));
            m_conf.read_config_line(s);
        }
        else if (a == "MV" || a == "fMV" ||  a == "GBW") {
            gdist_type = a;
        }
        else if (::isdigit(a[0])) {
            vector<string> pTnums = split(a, ",");
            for (vector<string>::iterator it = pTnums.begin(); it != pTnums.end(); it++) {
                m_pT.push_back(trim(*it, " \t"));
            }
        }
        else {
            // try opening as a file
            ifstream config;
            config.open(a.c_str());
            if (config.good()) {
                logger << "Reading config file " << a << endl;
                config >> m_conf;
                config.close();
            }
            else {
                m_conf.add("hardfactor_specifications", a);
            }
        }
    }

    if (!m_pT.empty()) {
        m_conf.erase("pT");
        for (vector<string>::iterator it = m_pT.begin(); it != m_pT.end(); it++) {
            m_conf.add("pT", *it);
        }
    }
    if (!gdist_type.empty()) {
        m_conf.set("gdist", gdist_type);
    }

    if (!m_conf.contains("hardfactor_specifications")) {
        m_conf.add("hardfactor_specifications", "lo");
        m_conf.add("hardfactor_specifications", "nlo");
    }
    pair<Configuration::iterator,Configuration::iterator> hf_bounds = m_conf.equal_range("hardfactor_specifications");
    for (Configuration::iterator it = hf_bounds.first; it != hf_bounds.second; it++) {
        m_hfspecs.push_back(it->first);
    }
}

ProgramConfiguration::~ProgramConfiguration() {
}

