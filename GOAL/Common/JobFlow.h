////////////////////////////////
/// usage : 1.	support running homogenous jobs on stream processors following CUDA's grid-stride loops paradigm.
/// 
/// note  : 1.	https://developer.nvidia.com/blog/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
////////////////////////////////

#ifndef CN_HUST_GOAL_JOB_FLOW_H
#define CN_HUST_GOAL_JOB_FLOW_H


#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <vector>


namespace goal {

namespace impl {

class JobFlowBase {
public:
    // index for thread/worker/job.
    // must be signed type!
    using ID = int;

    // `dim` is the total number of workers.
    // `idx` is the worker assigned to this job.
    using Job = std::function<void(ID dim, ID idx)>;


    // `init()` must be called before `run()` if `JobFlow` is instantiated by the trivial default constructor.
    virtual void init(ID workerNumber) = 0;
    virtual void init(ID workerNumber, Job kernel) {
        job = kernel;
        init(workerNumber);
    }

    // block the current thread until all jobs are done.
    virtual void run() = 0;
    virtual void run(Job kernel) {
        job = kernel;
        run();
    }

    ID workerNum() const { return sCast<ID>(homoWorkers.size()); }

protected:
    std::vector<std::thread> homoWorkers;
    Job job;
};


class SimpleJobFlow : public JobFlowBase {
public:
    SimpleJobFlow() {}
    SimpleJobFlow(ID workerNumber) : homoWorkers(workerNumber) {}

    using JobFlowBase::init;
    virtual void init(ID workerNumber) { homoWorkers.resize(workerNumber); }

    using JobFlowBase::run;
    virtual void run() {
        for (ID t = 0; t < workerNum(); ++t) {
            homoWorkers[t] = std::thread([&](ID index) { job(workerNum(), index); }, t);
        }
        for (auto t = homoWorkers.begin(); t != homoWorkers.end(); ++t) {
			if (t->joinable()) { t->join(); }
		}
    }

protected:
    Job job;
    std::vector<std::thread> homoWorkers;
};


#pragma warning(push)
#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
class JobFlowByFlag : public JobFlowBase {
public:
    JobFlowByFlag() {}
    JobFlowByFlag(ID workerNumber) { init(workerNumber); }

    ~JobFlowByFlag() {
        reset(-1);
        for (auto t = homoWorkers.begin(); t != homoWorkers.end(); ++t) {
            if (t->joinable()) { t->join(); }
        }
    }

    using JobFlowBase::init;
    virtual void init(ID workerNumber) {
        pendingWorkerNum = workerNumber;
        existPendingJobs.resize(workerNumber, 0);
        homoWorkers.reserve(workerNumber);
        for (ID t = 0; t < workerNumber; ++t) {
            homoWorkers.emplace_back([this](ID index) {
                for (;;) {
                    {
                        std::unique_lock<std::mutex> ul(m);
                        if (--pendingWorkerNum <= 0) { cvWait.notify_one(); } // all workers have done their jobs.
                        cvRun.wait(ul, [&]() { return existPendingJobs[index]; }); // wait for a new job for this worker.
                    }
                    if (pendingWorkerNum < 0) { return; } // all workers are fired.
                    job(workerNum(), index); // doing the job.
                    existPendingJobs[index] = 0; // job is done.
                }
            }, t);
        }
        wait();
    }

    using JobFlowBase::run;
    virtual void run() {
        reset(workerNum());
        wait();
    }

protected:
    void reset(ID num) {
        pendingWorkerNum = num;
        std::fill(existPendingJobs.begin(), existPendingJobs.end(), 1);
        cvRun.notify_all();
    }

    void wait() {
        std::unique_lock<std::mutex> ul(m);
        cvWait.wait(ul, [this]() { return pendingWorkerNum <= 0; });
    }

    ID pendingWorkerNum;
    std::vector<char> existPendingJobs; // do not use `std::vector<bool>`, it is not thread-safe even if you are accessing different index since they may share the same byte!

    std::mutex m;
    std::condition_variable cvRun;
    std::condition_variable cvWait;
};


class JobFlowByCounter : public JobFlowBase {
public:
    JobFlowByCounter() {}
    JobFlowByCounter(ID workerNumber) { init(workerNumber); }

    ~JobFlowByCounter() {
        reset(-1);
        for (auto t = homoWorkers.begin(); t != homoWorkers.end(); ++t) {
            if (t->joinable()) { t->join(); }
        }
    }

    using JobFlowBase::init;
    virtual void init(ID workerNumber) {
        pendingJobNum = 0;
        homoWorkers.reserve(workerNumber);
        for (ID t = 0; t < workerNumber; ++t) {
            homoWorkers.emplace_back([this](ID index) {
                for (;;) {
                    {
                        std::unique_lock<std::mutex> ul(m);
                        cvRun.wait(ul, [&]() { return pendingJobNum > 0; }); // wait for new jobs for all workers.
                        if (--pendingJobNum <= 0) { // picked a job.
                            ul.unlock();
                            cvRun.notify_all();
                        }
                    }
                    if (pendingWorkerNum < 0) { return; } // all workers are fired.
                    job(workerNum(), index); // doing the picked job.
                    {
                        std::unique_lock<std::mutex> ul(m);
                        cvRun.wait(ul, [&]() { return pendingJobNum <= 0; }); // wait until each worker has picked a job.
                        if (--pendingWorkerNum <= 0) { // all workers have done their jobs.
                            ul.unlock();
                            cvWait.notify_one();
                        }
                    }
                }
            }, t);
        }
    }

    using JobFlowBase::run;
    virtual void run() {
        reset(workerNum());
        wait();
    }

protected:
    void reset(ID num) {
        pendingWorkerNum = num;
        pendingJobNum = workerNum();
        cvRun.notify_all();
    }

    void wait() {
        std::unique_lock<std::mutex> ul(m);
        cvWait.wait(ul, [this]() { return pendingWorkerNum <= 0; });
    }

    ID pendingWorkerNum;
    ID pendingJobNum;

    std::mutex m;
    std::condition_variable cvRun;
    std::condition_variable cvWait;
};


class JobFlowByAccCounter : public JobFlowBase {
public:
    JobFlowByAccCounter() {}
    JobFlowByAccCounter(ID workerNumber) { init(workerNumber); }

    ~JobFlowByAccCounter() {
        reset(-1);
        for (auto t = homoWorkers.begin(); t != homoWorkers.end(); ++t) {
            if (t->joinable()) { t->join(); }
        }
    }

    using JobFlowBase::init;
    virtual void init(ID workerNumber) {
        pendingWorkerNum = workerNumber;
        accJobNum = 0;
        accPendingJobNums.resize(workerNumber, 0);
        homoWorkers.reserve(workerNumber);
        for (ID t = 0; t < workerNumber; ++t) {
            homoWorkers.emplace_back([this](ID index) {
                for (;;) {
                    {
                        std::unique_lock<std::mutex> ul(m);
                        if (--pendingWorkerNum <= 0) { cvWait.notify_one(); } // all workers have done their jobs.
                        cvRun.wait(ul, [&]() { return accPendingJobNums[index] < accJobNum; }); // wait for a new job for this worker.
                    }
                    if (pendingWorkerNum < 0) { return; } // all workers are fired.
                    job(workerNum(), index); // doing the job.
                    ++accPendingJobNums[index]; // job is done.
                }
            }, t);
        }
        wait();
    }

    using JobFlowBase::run;
    virtual void run() {
        reset(workerNum());
        wait();
    }

protected:
    void reset(ID num) {
        pendingWorkerNum = num;
        if (++accJobNum > (1 << 30)) {
            accJobNum = 1;
            std::fill(accPendingJobNums.begin(), accPendingJobNums.end(), 0);
        }
        cvRun.notify_all();
    }

    void wait() {
        std::unique_lock<std::mutex> ul(m);
        cvWait.wait(ul, [this]() { return pendingWorkerNum <= 0; });
    }

    ID pendingWorkerNum;
    ID accJobNum;
    std::vector<ID> accPendingJobNums; // do not use `std::vector<bool>`, it is not thread-safe even if you are accessing different index since they may share the same byte!

    std::mutex m;
    std::condition_variable cvRun;
    std::condition_variable cvWait;
};
#pragma warning(pop)

}


using JobFlow = impl::JobFlowByFlag;


namespace demo {

struct JobFlowDemo {
    static void demo0() {
        using ID = JobFlow::ID;

        ID threadNum = std::thread::hardware_concurrency();

        std::vector<int> a({ 1, 2, 3, 4, 5, 6 });

        JobFlow jobFlow;
        jobFlow.init(threadNum, [&](ID dim, ID idx) {
            ID size = sCast<ID>(a.size());
            for (ID i = idx; i < size; i += dim) { ++a[i]; }
        });
        jobFlow.run();
    }

    static void demo1() {
        using ID = JobFlow::ID;

        ID threadNum = std::thread::hardware_concurrency();

        std::vector<int> a({ 1, 2, 3, 4, 5, 6 });

        JobFlow jobFlow(threadNum);
        jobFlow.run([&](ID dim, ID idx) {
            ID size = sCast<ID>(a.size());
            for (ID i = idx; i < size; i += dim) { ++a[i]; }
        });
    }
};

}

}

#endif // CN_HUST_GOAL_JOB_FLOW_H
