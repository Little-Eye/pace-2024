////////////////////////////////
/// usage : 1.	definition of the data interface.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_OCM_PROBLEM_H
#define CN_HUST_GOAL_OCM_PROBLEM_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Optimization/Interface/EnvironmentBase.h"
#include "GOAL/Optimization/Interface/PerformanceBase.h"
#include "GOAL/Optimization/OCM/OcmTypedef.h"

#include <iostream>

namespace goal {
namespace OCM {

using Coord = Real;

#pragma warning(push)
#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
struct Input {
	static constexpr char const* TextInstanceFileExt = ".gr";

	// data for algorithm core.
	enum Simplify {
		NO,
		Neighbor
	};

	Simplify simplify = Simplify::Neighbor;

	Int nodeNum1, nodeNum2, arcNum;
	Vec<Vec<InNode>> inNodes;
	//Vec<NodeId> simplifyNodeId2GlobalNodeId;
	Vec<Vec<NodeId>> simplifyNodeIdSets;

	void simplifyByNeighbor();

	//Pack
	bool isPack = false;
	Int max_nodeNum2 = 20000;
	Vec<NodeIds> packNodeIdSets;

	void packByRandom();

	bool load(std::istream& ifs);
	bool load(const Str& path);
	bool load() { return load(std::cin); }

	void init() {
		switch (simplify)
		{
		case goal::OCM::Input::Neighbor:
			simplifyByNeighbor();
			break;
		}

		packByRandom();
	}
};

// solution vector.
struct Output {
	static constexpr char const* SolutionFileExt = ".sol";

	NodeIds order;
	CrossingNum cn;

	bool save(std::ostream& ofs, const Input& input) const;
	bool save(const Str& path, const Input& input) const;
	bool save(const Input& input) const { return save(std::cout, input); }
	bool load(const Str& path, const Input& input);


	Obj getObjValue() const { return cn; }
	Obj calcObjValue(const Str& inputPath, const Str& outputPath) const;
};
#pragma warning(pop)


enum ErrorFlag {
	Ok = 0x0,
	SizeError = 0x1,
	NumError = 0x2,

	// problem specific error.
	ProblemSpecificError = 0x4,
};
using ErrorFlags = int;

ErrorFlags check(const Input& input, const Output& output);
bool checkProblemSpecificConstraint(const Input& input, const Output& output);

ErrorFlags record(const Input& input, const Output& output, const Str& inputPath, const Str& outputPath, const EnvBase& env, const BenchmarkInfo& bi, const PerfBenchmark& perf, Record<Obj>& rec);

Records<Obj> loadBaseline(const Str& path);
void saveBaseline(const Str& path, const Records<Obj>& records);

}
}


#endif // CN_HUST_GOAL_OCM_PROBLEM_H
