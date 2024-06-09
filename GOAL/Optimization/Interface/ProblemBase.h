////////////////////////////////
/// usage : 1.	demo of the definition of the data interface.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_PROBLEM_H
#define CN_HUST_GOAL_COMMON_PROBLEM_H


#include "GOAL/Typedef.h"
#include "EnvironmentBase.h"
#include "PerformanceBase.h"


namespace goal {
namespace example {

using ObjectiveValue = double;

// benchmark instance.
// it is not necessary to inherit this class in order to implement an input.
struct Input {
    bool load(const Str& path) {}
};

// solution vector.
// it is not necessary to inherit this class in order to implement an output.
struct Output {
    ObjectiveValue getObjValue() const { return 0; }

    bool save(const Str& path, const Input& input) const;
    bool load(const Str& path, const Input& input);
};


enum ErrorFlag {
    Ok = 0x0,
    IoError = 0x1,
    FormatError = 0x2,
    InfeasibleError = 0x4,
};
using ErrorFlags = int;

ErrorFlags check(ObjectiveValue &obj, const Input &input, const Output &output);

ErrorFlags record(const Input &input, const Output &output, const EnvBase &env, const BenchmarkInfo &bi, const PerfBenchmark &perf);

}
}


#endif // CN_HUST_GOAL_COMMON_PROBLEM_H
