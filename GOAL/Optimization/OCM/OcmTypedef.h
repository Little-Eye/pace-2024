////////////////////////////////
/// usage : 1.	type aliases for simple types.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_OCM_TYPEDEF_H
#define CN_HUST_GOAL_COMMON_OCM_TYPEDEF_H

#include <utility>
#include <stack>
#include <math.h>
#include <bitset>
#include <queue>

#include "GOAL/Typedef.h"
#include "GOAL/Common/ConsecutiveIdSet.h"

using namespace std;

namespace goal {
	namespace OCM {
		using NodeId = Int;
		using Index = Int;
		using CrossingNum = LL;
		using Obj = CrossingNum;
		using Delta = Int;
		using Degree = Int;
		using NodeIds = Vec<NodeId>;

		struct InNode {
			NodeId inNodeId;
			Int num;

			InNode(NodeId _nodeId) {
				inNodeId = _nodeId;
				num = 1;
			}

			InNode(NodeId _nodeId, Int _num) {
				inNodeId = _nodeId;
				num = _num;
			}
		};

		struct Move {
			enum MoveFlag {
				UP,
				DOWN,
				NO
			};

			MoveFlag flag = MoveFlag::NO;
			Delta delta = 0;
			Index moveNodeId_index, beMovedNodeId_index;
		};
		using Moves = Vec<Move>;

		template<typename Id, typename Val>
		struct SortStruct {
			Id id;
			Val val;

			SortStruct() {
				id = -1;
				val = 0;
			}

			SortStruct(Id _id, Val _val) : id(_id), val(_val) {

			}

			bool operator<(const SortStruct& ss) const {
				return val < ss.val;
			}
		};

		struct Solution {
			CrossingNum cn;
			Vec<NodeId> order;
			Vec<Index> index;
			//Vec<Int> nodeTabuVec;
			Vec<Vec<Int>> tabuMat;

			void init(Int nodeNum) {
				cn = -1;
				order.resize(nodeNum);
				index.resize(nodeNum);
				tabuMat.resize(nodeNum);
				for (auto& tabuList : tabuMat) { tabuList.resize(nodeNum, -1); }
			}

			void clearTabuMat() {
				for (auto& tabuList : tabuMat) { fill(tabuList.begin(), tabuList.end(), -1); }
			}
		};
		using Solutions = Vec<Solution>;

		struct SCC {
			enum DynamicNodesBestImpLSFlag {
				dnBestImp,
				//dnPack,
				dnRestart
			};

			bool isSorted = false, isSmallSize = false;
			Int conNoImpEndIter = 0;
			Int nodeNum, moveLen, initNodeSize;

			//Tabu
			Int tabuMatFixedLen;
			Int tabuMatRandLen;

			//DynamicNodesBestImpLS
			Int dn_conBestImpIterLen;
			DynamicNodesBestImpLSFlag dynamicNodesBestImpLsFlag = DynamicNodesBestImpLSFlag::dnBestImp;

			//Solution
			Solution bestSol, searchSol;
			Int bestSol_size = 1;

			Vec<Vec<CrossingNum>> cm;

			Vec<NodeId> localNodeId2GlobalNodeId;

			void init(Int _nodeNum, Int baseMoveLen, Real tabuDividend, Int min_tabuMatFixedLen, Int min_tabuMatRandLen) {
				nodeNum = _nodeNum;
				initNodeSize = _nodeNum;
				moveLen = baseMoveLen;
				tabuMatFixedLen = max(min_tabuMatFixedLen, (Int)(_nodeNum / tabuDividend));
				//tabuMatRandLen = max(min_tabuMatRandLen, (Int)(_nodeNum / tabuDividend));
				tabuMatRandLen = 10;
				bestSol.init(nodeNum);
			}
		};

		struct TarjanData {
			Int dfncnt;
			Vec<Int> dfn, low;
			Vec<bool> inStack;
			std::stack<Int> sta;
			Int scccnt;
			Vec<Vec<Int>> scc;

			TarjanData(Int size) {
				dfncnt = 0;
				dfn.resize(size, -1), low.resize(size, -1);
				inStack.resize(size, false);
				scccnt = 0;
				scc.resize(size);
			}
		};
	}
}


#endif // CN_HUST_GOAL_COMMON_OCM_TYPEDEF_H
