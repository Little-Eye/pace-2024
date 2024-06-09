//#include <filesystem>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <atomic>
#include <chrono>
#include <utility>

#include "GOAL/Flag.h"
#include "GOAL/Common/Log.h"
#include "GOAL/Common/Timer.h"
#include "GOAL/Common/VarStr.h"
#include "GOAL/ThirdParty/Protobuf.h"
#include "GOAL/System/System.h"
#include "GOAL/System/File.h"
#include "GOAL/Optimization/OCM/Problem.h"
#include "GOAL/Optimization/OCM/Visualization.h"
//#include "GOAL/Optimization/OCM/ZstarBySzx.h"
//#include "GOAL/Optimization/OCM/DfsBySzx.h"
//#include "GOAL/Optimization/OCM/MipBySzx.h"
#include "GOAL/Optimization/OCM/OcmByKq.h"


// use multi-thread parallel benchmark instead of multi-process ones (e.g., job server like SGE, TORQUE, HTCondor).
#define SZX_MULTI_THREAD_BENCHMARK  1


using namespace std;
//using namespace std::experimental;
using namespace goal;
namespace ocm = OCM;
using Input = ocm::Input;
using Output = ocm::Output;
using Obj = ocm::Obj;
//using Alg = ocm::ZstarBySzx;
//using Alg = ocm::DfsBySzx;
using Alg = ocm::OcmByKq;

template<typename Algorithm = Alg>
bool run(const Str& inputPath, const Str& outputPath) {
	typename Algorithm::Environment env;
	env.msTimeout = 5 * 60 * 1000;
	env.maxThreadNum = 8;

	Input input;
	if (!input.load(inputPath)) { return false; }
	input.init();

	typename Algorithm::Configuration cfg;
	Algorithm alg;
	alg.init(input, env, cfg);

	Output output;
	alg.solve(output);

	output.save(outputPath, input);
	return true;
}

template<typename Algorithm = Alg>
bool solve() {
	typename Algorithm::Environment env;
	env.msTimeout = 5 * 60 * 1000;
	env.maxThreadNum = 1;
	env.randSeed = Random::generateSeed();

	Input input;
	if (!input.load()) { return false; }
	input.init();

	typename Algorithm::Configuration cfg;
	Algorithm alg;
	alg.init(input, env, cfg);

	Output output;
	alg.solve(output);

	output.save(input);
	return true;
}

template<typename Algorithm = Alg>
bool solve(Str inputPath, Str outputPath) {
	typename Algorithm::Environment env;
	env.msTimeout = 5 * 60 * 1000;
	env.maxThreadNum = 1;
	env.randSeed = Random::generateSeed();

	Input input;
	if (!input.load(inputPath)) { return false; }
	input.init();

	typename Algorithm::Configuration cfg;
	Algorithm alg;
	alg.init(input, env, cfg);

	Output output;
	alg.solve(output);

	output.save(outputPath, input);
	return true;
}


struct BenchmarkEnv : public Alg::Environment, public BenchmarkInfo {
	static Str InstanceDir() { return "OCM\\Instance\\heuristic-public\\"; }
	static Str SolutionDir() { return "OCM\\Solution\\heuristic-public\\"; }
	static Str VisualizDir() { return "OCM\\Visualization\\"; }
};

bool test(BenchmarkEnv& env, const Alg::Configuration& cfg, Record<Obj>& rec) {
	Stopwatch allSw;
	Stopwatch partSw;

	Log(Log::Level::Info) << "---- " << env.randSeed << endl;

	PerfBenchmark perf;
	Input input;
	if (!input.load(env.inputPath)) {
		Log(Log::Level::Error) << "failed to open instance " << env.inputPath << endl;
		return false;
	}                                                          partSw.printTime("load " + env.instanceName);
	input.init();

	Alg alg;
	alg.init(input, env, cfg, rec);                            partSw.printTime("init solver");

	Output output;
	perf.alg = alg.solve(output);                              partSw.printTime("run solver");

	perf.common.msTotalCpu = Timer::durationInMillisecond(allSw.lastTime, Timer::Clock::now());
	perf.common.peakMemory = os::peakMemoryUsage();

	#if !SZX_MULTI_THREAD_BENCHMARK
	this_thread::sleep_for(chrono::milliseconds(mt19937(random_device()())() % 8192));
	#endif // SZX_MULTI_THREAD_BENCHMARK
	env.cfgName = cfg.briefStr();
	output.save(env.outputPath, input);
	ocm::ErrorFlags err = ocm::record(input, output, env.inputPath, env.outputPath, env, env, perf, rec);			partSw.printTime("record results");

	/*if (obj <= rec.bestObj) {
		((env.outputPath += ("." + to_string(obj))) += ".") += Timer::getTightLocalTime();
		output.save(env.outputPath, input);
		if ((err == ocm::ErrorFlag::Ok) && math::updateMin(rec.bestObj, obj)) { rec.bestCpu = perf.alg.msConvergenceCpu; }
	}*/                                                          allSw.printTime(env.instanceName + " solved");

	return true;
}

void testAll(const Str& instanceListPath, const Str& baselinePath, const Str& logPath) {
	Records<Obj> records(ocm::loadBaseline(baselinePath));
	file::Lines instanceNames(file::readAllLines(instanceListPath, '#'));
	bool newRecords = false;

	BenchmarkEnv testEnv;
	testEnv.msTimeout = 5 * 60 * 1000; // TODO[szx][2]: set environment.
	testEnv.maxThreadNum = 8; // TODO[szx][2]: set environment.

	testEnv.id = "Y7000P2018H";
	testEnv.logPath = logPath;

	Alg::Configuration testCfg; // EXT[szx][3]: load configuration from file.

	//int benchmarkThreadNum = thread::hardware_concurrency()/4;
	int benchmarkThreadNum = 1;
	atomic<int> ai(0);
	Vec<thread> ts; ts.reserve(benchmarkThreadNum);
	auto run = [&]() {
		BenchmarkEnv env(testEnv); // EXT[szx][3]: load environment from file.
		Alg::Configuration cfg(testCfg);
		for (int i; (i = ai++) < sCast<int>(instanceNames.size());) {
			env.instanceName = instanceNames[i];
			env.inputPath = BenchmarkEnv::InstanceDir() + env.instanceName + Input::TextInstanceFileExt;
			env.outputPath = BenchmarkEnv::SolutionDir() + env.instanceName + Output::SolutionFileExt;
			env.startTime = Timer::getLocalTime();
			//env.randSeed = Random::generateSeed();
			env.randSeed = 1;

			//cfg.alg = Alg::Configuration::Algorithm::Trivial;
			if (test(env, cfg, records[env.instanceName])) { newRecords = true; }
		}
	};
	if (benchmarkThreadNum > 1) {
		for (int t = 0; t < benchmarkThreadNum; ++t) {
			ts.emplace_back(run);
			this_thread::sleep_for(chrono::seconds(4));
		}
		for (auto t = ts.begin(); t != ts.end(); ++t) { t->join(); }
	} else {
		run();
	}

	if (newRecords) { ocm::saveBaseline(baselinePath, records); }
}

void testOne(const Str& instanceName, const Str& baselinePath) {
	Records<Obj> records(ocm::loadBaseline(baselinePath));

	BenchmarkEnv env;
	//env.msTimeout = 8 * 3600 * 1000ll; // TODO[szx][2]: set environment.
	env.msTimeout = 0.00001 * 60 * 1000;
	//env.randSeed = Random::generateSeed(); // TODO[szx][2]: set environment.
	env.randSeed = 0;
	env.maxThreadNum = 1; // TODO[szx][2]: set environment.

	env.id = "HW.C";
	env.logPath = BenchmarkEnv::SolutionDir() + "0.log." + instanceName + ".csv";

	env.instanceName = instanceName;
	env.inputPath = BenchmarkEnv::InstanceDir() + env.instanceName + Input::TextInstanceFileExt;
	env.outputPath = BenchmarkEnv::SolutionDir() + env.instanceName + Output::SolutionFileExt;
	env.startTime = Timer::getLocalTime();

	Alg::Configuration cfg; // EXT[szx][3]: load configuration from file.
	test(env, cfg, records[env.instanceName]);
}

void checkSolution(const Str& inputPath, const Str& outputPath) {
	Input input;
	input.load(inputPath);
	Output output;
	output.load(outputPath, input);
	Obj obj;
	ostringstream errFlag;
	errFlag << "0x" << hex << setw(4) << setfill('0') << ocm::check(input, output);
	Log(Log::Level::Info) << "status=" << errFlag.str() << " obj=" << obj << endl;
}

void checkAllSolution(const Str& instanceListPath) {
	file::Lines instanceNames(file::readAllLines(instanceListPath, '#'));
	for (auto i = instanceNames.begin(); i != instanceNames.end(); ++i) {
		checkSolution(BenchmarkEnv::InstanceDir() + *i,
			BenchmarkEnv::SolutionDir() + *i + Output::SolutionFileExt);
	}
}

void drawSolution(const Str& inputPath, const Str& outputPath) {
	Input input;
	input.load(inputPath);
	Output output;
	output.load(outputPath, input);
	ocm::Visualizer::drawOutput(outputPath + ".html", input, output);
}


// syntax: `EXE {CmdOptionKey}{KeyValueDelim}VALUE`
// sample usage: `alg.exe -p=InstanceList.txt --LogPath=log.csv`.
// sample usage: `alg.exe "-i=case1.txt" -o="case1sln.txt"`.
struct CmdArg {
	enum CmdOptions {
		InstanceListPath,
		BaselinePath,
		LogPath,
		InstanceName,
		InstancePath,
		SolutionPath,

		Size,
	};

	static constexpr char KeyValueDelim = '=';

	const Map<Str, CmdOptions> CmdOptionKeys = {
		{ "--InstanceListPath", InstanceListPath },
		{ "-p", InstanceListPath },

		{ "--BaselinePath", BaselinePath },
		{ "-b", BaselinePath },

		{ "--LogPath", LogPath },
		{ "-l", LogPath },

		{ "--InstanceName", InstanceName },
		{ "-n", InstanceName },

		{ "--InstancePath", InstancePath },
		{ "-i", InstancePath },

		{ "--SolutionPath", SolutionPath },
		{ "-o", SolutionPath },
	};
	const Arr<Str, CmdOptions::Size> DefaultCmdOptions = {
		BenchmarkEnv::InstanceDir() + "0.InstanceList.txt",
		BenchmarkEnv::InstanceDir() + "0.Baseline.txt",
		BenchmarkEnv::SolutionDir() + "0.Log.csv",
		"",
		"",
		"",
	};

	Arr<Str, CmdOptions::Size> cmdOptions;

	static Str skipPrefix(const Str& str, const Str& prefix) { return str.substr(prefix.size() + 1); }
	static Str extractKey(const Str& str) { return str.substr(0, str.find(KeyValueDelim)); }

	bool parse(int argc, char* argv[]) {
		bool matchAnyKey = false;
		for (int i = 1; i < argc; ++i) {
			auto cmdOptionIndex = CmdOptionKeys.find(extractKey(argv[i])); // OPT[szx][9]: use prefix tree to parse.
			if (cmdOptionIndex == CmdOptionKeys.end()) { continue; }
			cmdOptions[cmdOptionIndex->second] = skipPrefix(argv[i], cmdOptionIndex->first);
			matchAnyKey = true;
		}

		for (int i = 0; i < CmdOptions::Size; ++i) {
			if (cmdOptions[i].empty()) { cmdOptions[i] = DefaultCmdOptions[i]; }
		}
		return matchAnyKey;
	}
};


int main(int argc, char* argv[]) {
	//filesystem::create_directories(BenchmarkEnv::SolutionDir);

	enum Flag {
		Test,
		Solve
	};

	Flag flag = Flag::Solve;
	CmdArg cmd;

	switch (flag)
	{
	case Test:
		if (cmd.parse(argc, argv) || (argc < 2)) {
			if (!cmd.cmdOptions[CmdArg::InstanceName].empty()) {
				testOne(cmd.cmdOptions[CmdArg::InstancePath], cmd.cmdOptions[CmdArg::SolutionPath]);
			}
			else {
				testAll(cmd.cmdOptions[CmdArg::InstanceListPath], cmd.cmdOptions[CmdArg::BaselinePath], cmd.cmdOptions[CmdArg::LogPath]);
			}
		}
		else {
			run<Alg>(argv[1], argv[2]);
		}
		break;
	case Solve:
		/*if (argc < 2) {
			solve();
		}
		else {
			solve(argv[1], argv[2]);
		}*/
		solve();
		break;
	default:
		break;
	}

	return 0;
}
