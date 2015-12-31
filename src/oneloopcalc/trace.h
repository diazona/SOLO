#pragma once

#include <bitset>

namespace trace_variable {
    enum {
        #define process(v) v,
        #include "../integration/ictx_var_list.inc"
        #undef process
        COUNT
    };

    // definition is in resultscalculator.cpp
    extern std::bitset<trace_variable::COUNT> trace_vars;
}

