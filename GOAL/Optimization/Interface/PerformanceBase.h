////////////////////////////////
/// usage : 1.	base classes for solver environments.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_PERFORMANCE_H
#define CN_HUST_GOAL_COMMON_PERFORMANCE_H


#include "GOAL/Typedef.h"
#include "GOAL/System/System.h"


namespace goal {

#pragma warning(push)
#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
// information to be recorded by solver for benchmark logging.
struct PerfWhiteBoxBase {
    Millisecond msConvergenceCpu; // duration in millisecond for finding the best solution.
    Str stopState; // can be iteration, gernation, tree depth, or node number.
};

// information to be recorded by benchmark code for benchmark logging.
struct PerfBlackBoxBase {
    Millisecond msTotalCpu; // duration in millisecond for finding the best solution.
    os::MemoryUsage peakMemory;
};
#pragma warning(pop)

struct PerfBenchmark {
    PerfWhiteBoxBase alg;
    PerfBlackBoxBase common;
};

}


#endif // CN_HUST_GOAL_COMMON_PERFORMANCE_H
