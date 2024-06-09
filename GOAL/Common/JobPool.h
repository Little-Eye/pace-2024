////////////////////////////////
/// usage : 1.	thread pool implemented by `std::atomic`.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_JOB_POOL_H
#define CN_HUST_GOAL_COMMON_JOB_POOL_H


#include <atomic>
#include <functional>
#include <thread>
#include <vector>


namespace goal {

class JobPool {
public:
	using ID = int;

	using IsJobTaken = std::function<bool()>;
	using Job = std::function<void(IsJobTaken isJobTaken)>;

	static void fightForJob(ID workerNum, Job userJob) {
		std::vector<std::thread> workers; workers.reserve(workerNum);
		std::atomic<ID> gt = 0; // global task index.
		for (ID w = 0; w < workerNum; ++w) {
			workers.emplace_back([&]() {
				ID t = 0; // local task index.
				userJob([&]() {
					ID ogt = gt; // old global task index.
					if (t < ogt) { ++t; return true; }
					return !gt.compare_exchange_strong(ogt, ++t, std::memory_order_relaxed);
				}); // https://en.cppreference.com/w/cpp/atomic/atomic/compare_exchange
			}); // https://en.cppreference.com/w/cpp/atomic/memory_order
		}
		for (ID w = 0; w < workerNum; ++w) { workers[w].join(); }
	}
};

namespace demo {

struct JobPoolDemo {
	static void demo0() {
		using ID = JobPool::ID;

		ID threadNum = std::thread::hardware_concurrency();

		std::vector<int> a({ 1, 2, 3, 4, 5, 6 });

		JobPool::fightForJob(threadNum, [&](JobPool::IsJobTaken isJobTaken) {
			for (auto i = a.begin(); i != a.end(); ++i) {
				if (!isJobTaken()) { ++(*i); }
			}
		});
	}
};

}

}


#endif // CN_HUST_GOAL_COMMON_JOB_POOL_H
