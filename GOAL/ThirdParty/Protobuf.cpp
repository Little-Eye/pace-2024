// for protobuf auto-linking.

#include "GOAL/Flag.h"


#if _PLUGIN_PROTOBUF && _CC_MS_VC

#if _DR_DEBUG
#pragma comment(lib, "libprotobufd")
#pragma comment(lib, "libprotobuf-lited")
#pragma comment(lib, "libprotocd")
#else // _DR_RELEASE
#pragma comment(lib, "libprotobuf")
#pragma comment(lib, "libprotobuf-lite")
#pragma comment(lib, "libprotoc")
#endif // _DR_DEBUG

#endif // _PLUGIN_PB
