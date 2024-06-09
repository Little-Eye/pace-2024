// for gurobi auto-linking.

#include "GOAL/Flag.h"


#if _PLUGIN_GUROBI

#include "Gurobi.h"

#include "GOAL/Common/Preprocessor.h"

#define GUROBI_VERSION  RESOLVED_STRINGIFY(RESOLVED_CONCAT(GRB_VERSION_MAJOR,GRB_VERSION_MINOR))

#if _DR_DEBUG
#if _LL_DYNAMIC
#define LINK_TYPE  "dd"
#else // _LL_STATIC
#define LINK_TYPE  "td"
#endif // _LL_DYNAMIC
#else // _DR_RELEASE
#if _LL_DYNAMIC
#define LINK_TYPE  "d"
#else // _LL_STATIC
#define LINK_TYPE  "t"
#endif // _LL_DYNAMIC
#endif // _DR_DEBUG

//#pragma message("[auto-linking] " "gurobi" GUROBI_VERSION)
//#pragma message("[auto-linking] " "gurobi_c++m" LINK_TYPE RESOLVED_STRINGIFY(_CC_VERSION))
#pragma comment(lib, "gurobi" GUROBI_VERSION)
#pragma comment(lib, "gurobi_c++m" LINK_TYPE RESOLVED_STRINGIFY(2017))


namespace goal {

thread_local GRBEnv MpSolverGurobi::globalEnv(true);
const MpSolverGurobi::Decision MpSolverGurobi::NoVar;

}

#endif // _PLUGIN_GUROBI
