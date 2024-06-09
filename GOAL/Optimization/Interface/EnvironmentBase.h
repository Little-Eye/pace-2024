////////////////////////////////
/// usage : 1.	base classes for solver environments.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_ENVIRONMENT_H
#define CN_HUST_GOAL_COMMON_ENVIRONMENT_H


#include <thread>
#include <limits>

#include "GOAL/Typedef.h"


namespace goal {

// the parameters which control initial program states, resource quota, and stop criteria.
struct EnvBase {
    virtual void calibrate() {
        int sysMaxThreadNum = std::thread::hardware_concurrency();
        if ((maxThreadNum <= 0) || (maxThreadNum > sysMaxThreadNum)) { maxThreadNum = sysMaxThreadNum; }
    }

    Str instanceName;
    Millisecond msTimeout; // run time limit in milliseconds.
    int randSeed; // seed of the random number generator.
    int maxThreadNum; // maximal available thread number.
};

struct EnvTreeSearchEnvBase {
    Iteration maxNodeNum; // maximal number of expanded tree nodes.
};
struct EnvTreeSearch : public EnvBase, public EnvTreeSearchEnvBase {};

struct EnvMetaheuristicsBase {
    Iteration maxIteration; // maximal number of performed neighborhood moves.
    Iteration maxStagnation; // maximal number of iterations without improving the local optima in the current trajectory.
    Iteration maxGeneration; // maximal number of crossovers or perturbations.
};
struct EnvMetaheuristics : public EnvBase, public EnvMetaheuristicsBase {};


// information to be set by benchmark code before solving for benchmark logging.
struct BenchmarkInfo {
    Str id; // can be platform id, hardware specification, or thread id.
    Str startTime; // time point `YYYY-MM-DD_hh:mm:ss`.
    Str cfgName;
    Str inputPath;
    Str outputPath;
    Str logPath;
};


// best result for each benchmark instance.
template<typename ObjectiveValue>
struct Record {
    static constexpr ObjectiveValue MaxObj = (std::numeric_limits<ObjectiveValue>::max)() / 2;
    static constexpr Millisecond MaxTime = (std::numeric_limits<Millisecond>::max)() / 2;

    ObjectiveValue refObj = MaxObj; // best known objective value in the literature.
    ObjectiveValue bestObj = MaxObj; // best known obective value including unpublished results.

    Millisecond refCpu = MaxTime; // CPU time to obtain the best known result in the literature.
    Millisecond bestCpu = MaxTime; // CPU time to obtain the best known result including unpublished results.

    ObjectiveValue bound = -MaxObj; // best known bound.
};
// best results for the problem.
template<typename ObjectiveValue>
using Records = Map<Str, Record<ObjectiveValue>>;

}


#endif // CN_HUST_GOAL_COMMON_ENVIRONMENT_H
