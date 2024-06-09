////////////////////////////////
/// usage : 1.	a simple hread pool without return value retrieval and argument passing.
/// 
/// note  : 1.	it is recommended to use `JobPool` instead.
////////////////////////////////

#ifndef HUST_HAFV_UTIL_THREAD_POOL_H
#define HUST_HAFV_UTIL_THREAD_POOL_H


#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <queue>
#include <functional>
#include <utility>

#include "GOAL/Typedef.h"


namespace goal {

namespace impl {

// [NoReturnValueRetrieval][NotExceptionSafe]
// [ManualStart][ManualPend][ManualStop]
class ThreadPoolBase {
public:
	enum State { Idle, Run, Pend, Stop };


	using Size = int;
	using Job = std::function<void(void)>;
	using Worker = std::thread;
	using Lock = std::unique_lock<std::mutex>;


	static Size getDefaultWorkerNum() { return Worker::hardware_concurrency(); }


	// you may want to start() once the pool is constructed in the derived classes.
	ThreadPoolBase(Size threadNum) : workerPool(threadNum) {}
	// you may want to pend() once the pool is destructed in the derived classes.
	virtual ~ThreadPoolBase() {}


	// make all workers ready to work.
	// you could restart after pend() or stop().
	virtual void start() {
		setState(State::Run);
		launchAll();
	}
	// terminate all workers after all pushed jobs are done.
	virtual void pend() {
		setState(State::Pend);
		waitAll(); // OPT[szx][0]: detach and move on instead of waiting?
	}
	// terminate all workers after all taken jobs are done.
	virtual void stop() {
		setState(State::Stop);
		waitAll();
	}

	virtual void push(Job&& newJob) = 0;

	virtual State getState() = 0;

protected:
	// the main loop of taking and executeing jobs for single worker thread.
	virtual void work() = 0;

	virtual void setState(State state) = 0;

	void launchAll() {
		for (auto worker = workerPool.begin(); worker != workerPool.end(); ++worker) {
			*worker = Worker([this]() { work(); });
		}
	}
	void waitAll() {
		for (auto worker = workerPool.begin(); worker != workerPool.end(); ++worker) {
			if (worker->joinable()) { worker->join(); }
		}
	}

	static bool isValidJob(const Job& job) { return sCast<bool>(job); }

	static void dummyJob() {}


	std::vector<Worker> workerPool;
};

// [NoReturnValueRetrieval][NotExceptionSafe]
// [ManualStart][ManualPend][ManualStop]
class ThreadPoolByQueue : public ThreadPoolBase {
public:
	using ThreadPoolBase::ThreadPoolBase;


	virtual void stop() override {
		setState(State::Stop);
		jobCv.notify_all(); // terminate sleeping workers.
		waitAll();
	}
	virtual void pend() override {
		setState(State::Pend);
		jobCv.notify_all(); // terminate sleeping workers.
		waitAll();
	}

	virtual void push(Job&& newJob) override {
		Lock jobLock(jobMutex);
		jobQueue.push(newJob);
		jobLock.unlock();

		jobCv.notify_one();
	}

	virtual State getState() override {
		Lock jobLock(jobMutex);
		return state;
	}

protected:
	virtual void work() override {
		for (;;) {
			Lock jobLock(jobMutex);

			while (jobQueue.empty()) {
				if (state != State::Run) { return; } // all pending jobs finished.
				jobCv.wait(jobLock);
			}
			if (state == State::Stop) { return; }

			Job newJob(std::move(jobQueue.front()));
			jobQueue.pop();

			jobLock.unlock();

			newJob();
		}
	}

	virtual void setState(State newState) override {
		Lock jobLock(jobMutex);
		state = newState;
	}


	State state;

	std::queue<Job> jobQueue;

	std::mutex jobMutex;
	std::condition_variable jobCv;
};


// [NoReturnValueRetrieval][NotExceptionSafe]
// [ManualStart][ManualPend][ManualStop]
class ThreadPoolBySingleSlot : public ThreadPoolBase {
public:
	using ThreadPoolBase::ThreadPoolBase;


	virtual void stop() override {
		state = State::Stop;
		jobCv.notify_all();
		waitAll();
	}
	virtual void pend() override {
		state = State::Pend;
		jobCv.notify_all();
		waitAll();
	}

	virtual void push(Job&& newJob) override {
		Lock workerLock(workerMutex);
		workerCv.wait(workerLock, [this]() { return isSlotEmpty(); }); // OPT[szx][0]: assume spurious wake up will never happen?
		nextJob = newJob; // make the new job available for taking.
		workerLock.unlock();

		jobCv.notify_one(); // declare that a new job is available (wake up a worker to take the new job).
	}

	virtual State getState() override { return state; }

protected:
	virtual void work() override {
		for (Job newJob; state != State::Stop; newJob()) { // do the job.
			Lock jobLock(jobMutex); // wait until a new job is available.

			while (isSlotEmpty()) {
				if (state != State::Run) { return; } // all pending jobs finished.
				jobCv.wait(jobLock);
			}
			if (state == State::Stop) { return; }

			newJob = std::move(nextJob); // take the new job.
			jobLock.unlock();

			workerCv.notify_one(); // declare that the job has been taken (wake up the dispatcher to enable a new job).
		}
	}

	virtual void setState(State newState) override { state = newState; }

	bool isSlotEmpty() const { return !isValidJob(nextJob); }


	std::atomic<State> state;

	Job nextJob;

	std::mutex jobMutex;
	std::condition_variable jobCv;
	std::mutex workerMutex;
	std::condition_variable workerCv;
};

}

// [NoReturnValueRetrieval][NotExceptionSafe]
// [AutoStart][AutoPend][ManualStop]
template<typename ThreadPoolImpl = impl::ThreadPoolByQueue>
class ThreadPool : public ThreadPoolImpl {
public:
	ThreadPool(int threadNum) : ThreadPoolImpl(threadNum) { ThreadPoolImpl::start(); }
	ThreadPool() : ThreadPool(ThreadPoolImpl::getDefaultWorkerNum()) {}
	virtual ~ThreadPool() { ThreadPoolImpl::pend(); }


	using ThreadPoolImpl::push;
	// avoid copying function objects. the const reference can be handled automatically.
	template<typename Functor>
	void push(Functor& newJob) { push(std::ref(newJob)); } // or use `push([&newJob]() { newJob(); });`.
};

}


#endif // HUST_HAFV_UTIL_THREAD_POOL_H
