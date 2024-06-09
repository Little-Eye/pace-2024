////////////////////////////////
/// usage : 1.	utilities for multi-threading.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_SYSTEM_THREAD_H
#define CN_HUST_GOAL_SYSTEM_THREAD_H


#include <atomic>
#include <thread>


namespace goal {

struct Thread {
	static std::atomic<int> gid; // global thread id.
	static thread_local int id; // thread id.

	static int num() { return std::thread::hardware_concurrency(); }
};

}


#endif // CN_HUST_GOAL_SYSTEM_THREAD_H
