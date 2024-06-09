// for zlib auto-linking.

#include "GOAL/Flag.h"


#if _PLUGIN_ZLIB && _CC_MS_VC

#if _DR_DEBUG
#pragma comment(lib, "zlibstaticd")
//#pragma comment(lib, "zlibd")
#else // _DR_RELEASE
#pragma comment(lib, "zlibstatic")
//#pragma comment(lib, "zlib")
#endif // _DR_DEBUG

#endif // _PLUGIN_PB
