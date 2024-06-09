////////////////////////////////
/// usage : 1.	sample code for solvers.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_SOLVER_BASE_H
#define CN_HUST_GOAL_COMMON_SOLVER_BASE_H


#include "ProblemBase.h"
#include "ConfigurationBase.h"
#include "PerformanceBase.h"


namespace goal {
namespace example {

// sample code for a basic solver.
// the basic interfaces must be implemented, and additional public methods are allowed.
// it is not necessary to inherit this class in order to implement a solver.
class SolverBase {
public:
    struct Environment : public EnvBase {

    };

    struct Configuration : public CfgBase {

    };

    struct Performance : public PerfWhiteBoxBase {

    };


	bool init(const Input& input, const Environment& env, 
        const Configuration& cfg = Configuration(),
        const Record<ObjectiveValue>& rec = Record<ObjectiveValue>()) { return false; }
	bool init(Input&& input, const Environment& env,
        const Configuration& cfg = Configuration(),
        const Record<ObjectiveValue>& rec = Record<ObjectiveValue>()) { return false; }
    Performance solve(Output &output) { return Performance(); }
};

// sample code for an incremental solver.
// the basic interfaces must be implemented, and additional public methods are allowed.
// it is not necessary to inherit this class in order to implement a solver.
class IncrementalSolverBase : public SolverBase {
public:
    template<typename T>
    bool update(const T &inputModifier) { return false; }
};

}
}


#endif // CN_HUST_GOAL_COMMON_SOLVER_BASE_H
