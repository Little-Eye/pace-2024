////////////////////////////////
/// usage : 1.	system signal handling.
/// 
/// note  : 1.	https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/signal?view=msvc-170
///         2.	https://docs.microsoft.com/en-us/previous-versions/ms811896(v=msdn.10)?redirectedfrom=MSDN#signals-and-signal-handling
///	        3.	https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/raise?view=msvc-170
////////////////////////////////

#ifndef CN_HUST_GOAL_SYSTEM_SIGNAL_H
#define CN_HUST_GOAL_SYSTEM_SIGNAL_H


#include <csignal>


namespace goal {
namespace os {

using SignalHandler = void (*)(int);


enum Signal {
	Default = 0, // SIG_DFL.
	Interrupt = SIGINT,
	IllegalInstruction = SIGILL,
	FloatingPointException = SIGFPE,
	SegmentFault = SIGSEGV,
	Termination = SIGTERM,
	Break = SIGBREAK,
	Abort = SIGABRT,
};


inline void setSignalHandler(SignalHandler signalHandler) {
	signal(SIGINT, signalHandler);
	signal(SIGBREAK, signalHandler);
	signal(SIGABRT, signalHandler);
	signal(SIGTERM, signalHandler);
}

inline void resetSignalHandler() {
	signal(SIGINT, SIG_DFL);
	signal(SIGBREAK, SIG_DFL);
	signal(SIGABRT, SIG_DFL);
	signal(SIGTERM, SIG_DFL);
}

inline void rasieSignal(int signal) { raise(signal); }

}
}


#endif // CN_HUST_GOAL_SYSTEM_SIGNAL_H
