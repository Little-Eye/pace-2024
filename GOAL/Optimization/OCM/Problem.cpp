#include "GOAL/Optimization/OCM/Problem.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mutex>
#include <cmath>

#include "GOAL/Common/Log.h"
#include "GOAL/Common/Container.h"
#include "GOAL/Common/Random.h"
#include "GOAL/Common/Math.h"
#include "GOAL/System/File.h"


using namespace std;


namespace goal {
namespace OCM {

#pragma region Input

void Input::simplifyByNeighbor() {
	Int local_nodeNum = 0;
	NodeIds simplifyNodeId2GlobalNodeId(nodeNum2, -1);
	simplifyNodeIdSets.resize(nodeNum2);

	for (NodeId global_unsimpNodeId = 0; global_unsimpNodeId < nodeNum2; ++global_unsimpNodeId) {
		bool isEqual = false;
		NodeId local_matchNodeId = -1;

		for (NodeId local_simpNodeId = 0; local_simpNodeId < local_nodeNum; ++local_simpNodeId) {
			NodeId global_simpNodeId = simplifyNodeId2GlobalNodeId[local_simpNodeId];
			if (inNodes[global_unsimpNodeId].size() != inNodes[global_simpNodeId].size()) { continue; }

			bool flag = false;
			Int matchNum = 0;
			for (auto neiNode1 : inNodes[global_unsimpNodeId]) {
				for (auto neiNode2 : inNodes[global_simpNodeId]) {
					if (neiNode1.inNodeId == neiNode2.inNodeId) {
						flag = true;
						matchNum++;
						break;
					}
				}
				if (!flag) { break; }
				flag = false;
			}

			if (matchNum == inNodes[global_unsimpNodeId].size()) {
				isEqual = true;
				local_matchNodeId = local_simpNodeId;
				break;
			}
		}

		if (!isEqual) {
			NodeId local_unsimpNodeId = local_nodeNum++;

			simplifyNodeIdSets[local_unsimpNodeId].emplace_back(global_unsimpNodeId);
			simplifyNodeId2GlobalNodeId[local_unsimpNodeId] = global_unsimpNodeId;
		}
		else {
			simplifyNodeIdSets[local_matchNodeId].emplace_back(global_unsimpNodeId);
		}
	}
	for (NodeId localNodeId = 0; localNodeId < local_nodeNum; ++localNodeId) {
		NodeId globalNodeId = simplifyNodeId2GlobalNodeId[localNodeId];
		inNodes[localNodeId] = move(inNodes[globalNodeId]);
		for (auto& inNode : inNodes[localNodeId]) { inNode.num = simplifyNodeIdSets[localNodeId].size(); }
	}

	nodeNum2 = local_nodeNum;
}

void Input::packByRandom() {
	if (nodeNum2 <= max_nodeNum2) { return; }

	isPack = true;
	packNodeIdSets.resize(max_nodeNum2);
	Int packNode_size = ceil((Real)nodeNum2 / max_nodeNum2);

	Int packNodeId_size_nodeNum = nodeNum2 - (packNode_size - 1) * max_nodeNum2;
	NodeId unpackNodeId = 0, packNodeId = 0;
	for (; packNodeId < packNodeId_size_nodeNum; ++packNodeId) {
		Int pPackNode_size = packNode_size;
		while (pPackNode_size-- && unpackNodeId < nodeNum2) {
			packNodeIdSets[packNodeId].emplace_back(unpackNodeId++);
		}
	}
	for (; packNodeId < max_nodeNum2; ++packNodeId) {
		Int pPackNode_size = packNode_size - 1;
		while (pPackNode_size-- && unpackNodeId < nodeNum2) {
			packNodeIdSets[packNodeId].emplace_back(unpackNodeId++);
		}
	}

	Vec<Vec<InNode>> pInNodes = move(inNodes);
	inNodes.resize(max_nodeNum2);
	NodeIds neiNodeId_buckets(nodeNum1);
	for (NodeId packNodeId = 0; packNodeId < max_nodeNum2; ++packNodeId) {
		fill(neiNodeId_buckets.begin(), neiNodeId_buckets.end(), 0);
		for (auto unpackNodeId : packNodeIdSets[packNodeId]) {
			for (auto& inNode : pInNodes[unpackNodeId]) { neiNodeId_buckets[inNode.inNodeId] += inNode.num; }
		}
		for (NodeId neiNodeId = 0; neiNodeId < nodeNum1; ++neiNodeId) {
			if (neiNodeId_buckets[neiNodeId] > 0) { inNodes[packNodeId].push_back({ neiNodeId, neiNodeId_buckets[neiNodeId] }); }
		}
	}

	nodeNum2 = max_nodeNum2;
}

bool Input::load(istream& ifs) {
	Str inputStr = "";
	Str p = "", ocr = "";
	getline(ifs, inputStr);
	istringstream iss(inputStr);

	iss >> p;
	while (p[0] == 'c') {		//check comment line
		getline(ifs, inputStr);
		iss.clear();
		iss.str(inputStr);
		iss >> p;
	}
	iss >> ocr >> nodeNum1 >> nodeNum2 >> arcNum;
	inNodes.resize(nodeNum2);

	Int nodeId1, nodeId2;
	for (Int arcId = 0; arcId < arcNum; ++arcId) {
		getline(ifs, inputStr);
		iss.clear();
		iss.str(inputStr);
		iss >> nodeId1 >> nodeId2;
		nodeId1 = nodeId1 - 1, nodeId2 = nodeId2 - nodeNum1 - 1;		//both start from 0
		inNodes[nodeId2].push_back({ nodeId1 });
	}

	return true;
}

bool Input::load(const Str& inputPath) {
	ifstream ifs(inputPath);
	if (!ifs.is_open()) { return false; }

	return load(ifs);
}
#pragma endregion Input


ErrorFlags check(const Input& input, const Output& output) {
	ErrorFlags error = ErrorFlag::Ok;

	return error;
}

#pragma region Check
bool checkProblemSpecificConstraint(const Input& input, const Output& output) {
	return true;
}
#pragma endregion Check

#pragma warning(push)
#pragma warning(disable: 26115) // Warning C26115 Failing to release lock.
ErrorFlags record(const Input& input, const Output& output, const Str& inputPath, const Str& outputPath, const EnvBase& env, const BenchmarkInfo& bi, const PerfBenchmark& perf, Record<Obj>& rec) {
	ErrorFlags err = check(input, output);

	// format the log.
	ostringstream log;
	ostringstream errFlag;
	errFlag << "0x" << hex << setw(4) << setfill('0') << err;
	log << bi.startTime << "," << bi.id << "," << env.instanceName << "," << errFlag.str() << "," <<
		0 << "," << output.calcObjValue(inputPath, outputPath) << ","
		<< perf.alg.msConvergenceCpu << "," << perf.common.peakMemory.physicalMemory << "," << perf.common.peakMemory.virtualMemory << ","
		<< env.randSeed << "," << bi.cfgName << "," << perf.alg.stopState << ",";
	
	// EXT[szx][9]: record solution vector.
	log << endl;

	// append the log.
	static mutex logFileMutex;
	lock_guard<mutex> logFileGuard(logFileMutex);

	ofstream logFile(bi.logPath, ios::app);
	logFile.seekp(0, ios::end);
	if (logFile.tellp() <= 0) {
		logFile << "Time,ID,Instance,Feasible,ObjErr,Obj,Duration,PhysMem,VirtMem,RandSeed,Config,Stop,Solution" << endl;
	}
	logFile << log.str();
	logFile.close();
	return err;
}
#pragma warning(pop)

bool Output::save(ostream& ofs, const Input& input) const {
	for (auto nodeId : order) {
		ofs << nodeId + input.nodeNum1 + 1 << endl;
	}

	return true;
}

bool Output::save(const Str& outputPath, const Input& input) const {
	ofstream ofs(outputPath);
	if (!ofs.is_open()) { return false; }

	bool ret = save(ofs, input);

	ofs.close();
	return ret;
}

bool Output::load(const Str& outputPath, const Input& input) {
	Str charBuf(file::readAllText(outputPath));
	if (charBuf.empty()) { charBuf = file::readAllText(outputPath); }
	if (charBuf.empty()) {
		Log(Log::Level::Fatal) << "failed to open solution file: " << outputPath << endl;
		return false;
	}
	istringstream ifs(charBuf);

	return true;
}

Obj Output::calcObjValue(const Str& inputPath, const Str& outputPath) const {
	ifstream ifs(inputPath);
	Str inputStr;
	Str p, ocr;
	Int nodeNum1, nodeNum2, arcNum;
	Vec<NodeIds> in_nodeIds;
	getline(ifs, inputStr);
	istringstream iss(inputStr);

	iss >> p;
	while (p[0] == 'c') {		//check comment line
		getline(ifs, inputStr);
		iss.clear();
		iss.str(inputStr);
		iss >> p;
	}
	iss >> ocr >> nodeNum1 >> nodeNum2 >> arcNum;
	in_nodeIds.resize(nodeNum2);

	Int nodeId1, nodeId2;
	for (Int arcId = 0; arcId < arcNum; ++arcId) {
		getline(ifs, inputStr);
		iss.clear();
		iss.str(inputStr);
		iss >> nodeId1 >> nodeId2;
		nodeId1 = nodeId1 - 1, nodeId2 = nodeId2 - nodeNum1 - 1;		//both start from 0
		in_nodeIds[nodeId2].emplace_back(nodeId1);
	}

	Obj cn = 0;
	NodeIds order(nodeNum2);
	ifstream ofs(outputPath);
	for (Int i = 0; i < nodeNum2; ++i) { ofs >> order[i]; }
	for (auto nodeId1 = order.begin(); nodeId1 != order.end(); ++nodeId1) {
		for (auto nodeId2 = nodeId1 + 1; nodeId2 != order.end(); ++nodeId2) {
			NodeId t_nodeId1 = *nodeId1 - nodeNum1 - 1, t_nodeId2 = *nodeId2 - nodeNum1 - 1;
			for (auto nei_nodeId1 : in_nodeIds[t_nodeId1]) {
				for (auto nei_nodeId2 : in_nodeIds[t_nodeId2]) {
					if (nei_nodeId1 > nei_nodeId2) { cn++; }
				}
			}
		}
	}

	return cn;
}

Records<Obj> loadBaseline(const Str& path) {
	Records<Obj> records;
	ifstream ifs(path);
	Str s;
	for (getline(ifs, s); getline(ifs, s);) { // skip the header.
		istringstream iss(s);
		iss >> s;
		iss >> records[s].bestObj >> records[s].refObj >> records[s].bestCpu >> records[s].refCpu >> records[s].bound;
	}
	return records;
}
void saveBaseline(const Str& path, const Records<Obj>& records) {
	static constexpr char Delimiter = '\t';

	ofstream ofs(path);
	ofs << "Instance" << Delimiter << "BestObj" << Delimiter << "RefObj" << Delimiter << "BestCpu" << Delimiter << "RefCpu" << Delimiter << "Bound" << endl;
	for (auto r = records.begin(); r != records.end(); ++r) {
		ofs << r->first << Delimiter << r->second.bestObj << Delimiter << r->second.refObj
			<< Delimiter << r->second.bestCpu << Delimiter << r->second.refCpu << Delimiter << r->second.bound << endl;
	}
}

}
}
