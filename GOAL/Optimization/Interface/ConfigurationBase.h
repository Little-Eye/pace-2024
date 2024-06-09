////////////////////////////////
/// usage : 1.	base classes for solver configurations.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_CONFIGURATION_H
#define CN_HUST_GOAL_COMMON_CONFIGURATION_H


#include "GOAL/Typedef.h"


namespace goal {

// the parameters which control the performance of the algorithms.
struct CfgBase {
    virtual Str briefStr() const { return "NULL"; }
};

struct CfgPopulationBasedMetaheuristicsBase {
    Iteration maxPopulation; // maximal number of individuals in the population.
};
struct CfgPopulationBasedMetaheuristics : public CfgBase, public CfgPopulationBasedMetaheuristicsBase {};

}


#endif // CN_HUST_GOAL_COMMON_CONFIGURATION_H
