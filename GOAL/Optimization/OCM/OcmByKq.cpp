#include "GOAL/Optimization/OCM/OcmByKq.h"

#include <sstream>
#include <signal.h>

using namespace std;

namespace goal {
	namespace OCM {
		volatile sig_atomic_t isClosedBySig = 0;

		#pragma region Signal
		void sigtermHandler(Int signum) {
			isClosedBySig = 1;
		}
		#pragma endregion Signal

		Str OcmByKq::Configuration::briefStr() const {
			std::ostringstream oss;
			switch (initAlg)
			{
			case Configuration::InitAlgorithm::IndegreeSort:
				oss << "IndegreeSort ";
				break;
			case Configuration::InitAlgorithm::CrossingSort :
				oss << "CrossingSort";
				break;
			case Configuration::InitAlgorithm::UpAndDownCrossingSort :
				oss << "UpAndDownCrossingSort";
				break;
			case Configuration::InitAlgorithm::QuickSort :
				oss << "QuickSort";
				break;
			case Configuration::InitAlgorithm::MergeSort :
				oss << "MergeSort";
				break;
			}
			switch (alg) {
			case Algorithm::DynamicNodesBestImpLocalSearch :
				oss << "DynamicNodesBestImp";
				break;
			}
			return oss.str();
		}

		#pragma region basic
		bool OcmByKq::init(const Input& input, const Environment& env,
			const Configuration& cfg,
			const Record<CrossingNum>& rec) {
			this->pInput = &input;
			this->env = env;
			this->cfg = cfg;
			rand = Random(env.randSeed);
			timer = Timer(env.msTimeout);
			isClosedBySig = false;

			nodeNum = input.nodeNum2;
			initCM();

			initSCCs();

			initSolution();

			return true;
		}

		OcmByKq::Performance OcmByKq::solve(Output& output) {
			Performance retPer;

			switch (cfg.alg) {
			case Configuration::Algorithm::DynamicNodesBestImpLocalSearch :
				retPer = solve_DynamicNodes_BestImp_LocalSearch(output);
				break;
			}

			//visualize
			/*Str inputFileName = "OCM\\Instance\\heuristic-public\\Pic\\" + env.instanceName, outputFileName = "OCM\\Solution\\heuristic-public\\Pic\\" + env.instanceName;
			Visualizer visualizer;
			visualizer.drawInput(inputFileName, *pInput);
			visualizer.drawOutput(outputFileName, *pInput, output);*/

			return retPer;
		}

		bool OcmByKq::isBreakOut() {
			return isClosedBySig || timer.isTimeOut();
		}

		void OcmByKq::saveOutput(Vec<SCC>& sccs, Output& output) {
			output.order.reserve(nodeNum);

			if (pInput->isPack) {
				for (auto& scc : sccs) {
					for (auto localNodeId : scc.bestSol.order) {
						NodeId pack_globalNodeId = scc.localNodeId2GlobalNodeId[localNodeId];
						for (auto unpack_globalNodeId : (pInput->packNodeIdSets)[pack_globalNodeId]) {
							output.order.emplace_back(unpack_globalNodeId);
						}
					}
				}
			}
			else {
				for (auto& scc : sccs) {
					for (auto localNodeId : scc.bestSol.order) {
						NodeId unpack_globalNodeId = scc.localNodeId2GlobalNodeId[localNodeId];
						output.order.emplace_back(unpack_globalNodeId);
					}
				}
			}

			switch (pInput->simplify)
			{
			case Input::Simplify::Neighbor :
				NodeIds pOrder = move(output.order);
				for (auto simp_globalNodeId : pOrder) {
					for (auto unsimp_globalNodeId : (pInput->simplifyNodeIdSets)[simp_globalNodeId]) {
						output.order.emplace_back(unsimp_globalNodeId);
					}
				}
			}

			/*switch (pInput->simplify)
			{
			case Input::Simplify::NO:
				for (auto& scc : sccs) {
					for (auto localNodeId : scc.bestSol.order) {
						output.order.emplace_back(scc.localNodeId2GlobalNodeId[localNodeId]);
					}
				}
				break;
			case Input::Simplify::Neighbor:
				for (auto& scc : sccs) {
					for (auto localNodeId : scc.bestSol.order) {
						NodeId globalPreNodeId = scc.localNodeId2GlobalNodeId[localNodeId];
						for (auto globalNodeId : (pInput->simplifyNodeIdSets)[globalPreNodeId]) {
							output.order.emplace_back(globalNodeId);
						}
					}
				}
				break;
			}*/
		}
		#pragma endregion basic
		
		#pragma region initial

		void OcmByKq::initIndegreeSort(SCC& scc) {
			Vec<Degree> indegree(scc.nodeNum, 0);
			Vec<SortStruct<NodeId, Int>> sortVec(scc.nodeNum);
			for (NodeId nodeId1 = 0; nodeId1 < scc.nodeNum; ++nodeId1) {
				for (NodeId nodeId2 = nodeId1 + 1; nodeId2 < scc.nodeNum; ++nodeId2) {
					if (scc.cm[nodeId1][nodeId2] < scc.cm[nodeId2][nodeId1]) { indegree[nodeId2]++; }		//nodeId1 wants to precede nodeId2, so indegree[nodeId2]++
					else if (scc.cm[nodeId1][nodeId2] > scc.cm[nodeId2][nodeId1]) { indegree[nodeId1]++; }
				}
			}
			for (Int i = 0; i < scc.nodeNum; ++i) {
				sortVec[i].id = i;
				sortVec[i].val = indegree[i];
			}

			sort(sortVec.begin(), sortVec.end());
			
			for (Int i = 0; i < scc.nodeNum; ++i) {
				scc.bestSol.order[i] = sortVec[i].id;
				scc.bestSol.index[sortVec[i].id] = i;
			}

			calculateCN(scc.cm, scc.bestSol);
		}

		void OcmByKq::initCrossingSort(SCC& scc) {
			ConsecutiveIdSet valid_nodeIds(scc.nodeNum);
			for (NodeId nodeId = 0; nodeId < scc.nodeNum; ++nodeId) { valid_nodeIds.insert(nodeId); }
			Vec<CrossingNum> crossings(scc.nodeNum, 0);		//crossings[nodeId] is the crossing when nodeId preceeds the remaining nodes
			for (NodeId nodeId1 = 0; nodeId1 < scc.nodeNum; ++nodeId1) {
				for (NodeId nodeId2 = nodeId1 + 1; nodeId2 < scc.nodeNum; ++nodeId2) {
					crossings[nodeId1] += scc.cm[nodeId1][nodeId2];
					crossings[nodeId2] += scc.cm[nodeId2][nodeId1];
				}
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				CrossingNum minCrossing;
				NodeId minCrossing_nodeId = -1;
				Int minCrossing_size = 0;
				for (Index valid_nodeId_index = 0; valid_nodeId_index < valid_nodeIds.size(); ++valid_nodeId_index) {
					NodeId valid_nodeId = valid_nodeIds.itemAt(valid_nodeId_index);
					if (minCrossing_size == 0 || minCrossing > crossings[valid_nodeId]) {
						minCrossing = crossings[valid_nodeId];
						minCrossing_nodeId = valid_nodeId;
						minCrossing_size = 1;
					}
					else if (minCrossing == crossings[valid_nodeId] && rand.isPicked(1, ++minCrossing_size)) {
						minCrossing_nodeId = valid_nodeId;
					}
				}

				scc.bestSol.order[index] = minCrossing_nodeId;
				scc.bestSol.index[minCrossing_nodeId] = index;
				valid_nodeIds.eraseItem(minCrossing_nodeId);

				//update crossings of the remaining nodes
				for (Index valid_nodeId_index = 0; valid_nodeId_index < valid_nodeIds.size(); ++valid_nodeId_index) {
					NodeId valid_nodeId = valid_nodeIds.itemAt(valid_nodeId_index);
					crossings[valid_nodeId] -= scc.cm[valid_nodeId][minCrossing_nodeId];
				}
			}

			calculateCN(scc.cm, scc.bestSol);
		}

		void OcmByKq::initUpCrossingAndDownCrossingSort(SCC& scc) {
			ConsecutiveIdSet valid_nodeIds(scc.nodeNum);
			for (NodeId nodeId = 0; nodeId < scc.nodeNum; ++nodeId) { valid_nodeIds.insert(nodeId); }
			Vec<CrossingNum> upCrossings(scc.nodeNum, 0);		//upCrossings[nodeId] is the crossing when nodeId preceeds the remaining nodes
			Vec<CrossingNum> downCrossings(scc.nodeNum, 0);		//downCrossings[nodeId] is the crossing when nodeId follow the remaining nodes
			for (NodeId nodeId1 = 0; nodeId1 < scc.nodeNum; ++nodeId1) {
				for (NodeId nodeId2 = nodeId1 + 1; nodeId2 < scc.nodeNum; ++nodeId2) {
					upCrossings[nodeId1] += scc.cm[nodeId1][nodeId2];
					upCrossings[nodeId2] += scc.cm[nodeId2][nodeId1];

					downCrossings[nodeId1] += scc.cm[nodeId2][nodeId1];
					downCrossings[nodeId2] += scc.cm[nodeId1][nodeId2];
				}
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				CrossingNum minVal;
				NodeId minVal_nodeId = -1;
				Int minVal_size = 0;
				for (Index valid_nodeId_index = 0; valid_nodeId_index < valid_nodeIds.size(); ++valid_nodeId_index) {
					NodeId valid_nodeId = valid_nodeIds.itemAt(valid_nodeId_index);
					CrossingNum val = upCrossings[valid_nodeId] - downCrossings[valid_nodeId];
					if (minVal_size == 0 || minVal > val) {
						minVal = val;
						minVal_nodeId = valid_nodeId;
						minVal_size = 1;
					}
					else if (minVal == val && rand.isPicked(1, ++minVal_size)) {
						minVal_nodeId = valid_nodeId;
					}
				}

				scc.bestSol.order[index] = minVal_nodeId;
				scc.bestSol.index[minVal_nodeId] = index;
				valid_nodeIds.eraseItem(minVal_nodeId);

				//update crossings of the remaining nodes
				for (Index valid_nodeId_index = 0; valid_nodeId_index < valid_nodeIds.size(); ++valid_nodeId_index) {
					NodeId valid_nodeId = valid_nodeIds.itemAt(valid_nodeId_index);
					upCrossings[valid_nodeId] -= scc.cm[valid_nodeId][minVal_nodeId];
					downCrossings[valid_nodeId] -= scc.cm[minVal_nodeId][valid_nodeId];
				}
			}
		}

		void OcmByKq::quickSortBlock(NodeIds& initNodeIds, Vec<Index>& initIndexs, SCC& scc) {
			Int initSize = initNodeIds.size();
			if (initSize == 0) { return; }
			else if (initSize == 1) {
				NodeId initNodeId = initNodeIds[0];
				Index initIndex = initIndexs[0];
				scc.bestSol.order[initIndex] = initNodeId;
				scc.bestSol.index[initNodeId] = initIndex;
				return;
			}
			else if (initSize == 2) {
				NodeId initNodeId1 = initNodeIds[0], initNodeId2 = initNodeIds[1];
				Index initUpIndex = initIndexs[0], initDownIndex = initIndexs[1];

				if (scc.cm[initNodeId1][initNodeId2] <= scc.cm[initNodeId2][initNodeId1]) {
					scc.bestSol.order[initUpIndex] = initNodeId1;
					scc.bestSol.index[initNodeId1] = initUpIndex;
					scc.bestSol.order[initDownIndex] = initNodeId2;
					scc.bestSol.index[initNodeId2] = initDownIndex;
				}
				else {
					scc.bestSol.order[initUpIndex] = initNodeId2;
					scc.bestSol.index[initNodeId2] = initUpIndex;
					scc.bestSol.order[initDownIndex] = initNodeId1;
					scc.bestSol.index[initNodeId1] = initDownIndex;
				}

				return;
			}

			Vec<Int> indegrees_nodeIdIndex(initSize, 0);
			for (Index initNodeId1_index = 0; initNodeId1_index < initSize; ++initNodeId1_index) {
				NodeId initNodeId1 = initNodeIds[initNodeId1_index];
				for (Index initNodeId2_index = initNodeId1_index + 1; initNodeId2_index < initSize; ++initNodeId2_index) {
					NodeId initNodeId2 = initNodeIds[initNodeId2_index];

					if (scc.cm[initNodeId1][initNodeId2] < scc.cm[initNodeId2][initNodeId1]) { indegrees_nodeIdIndex[initNodeId2_index]++; }
					else if (scc.cm[initNodeId1][initNodeId2] > scc.cm[initNodeId2][initNodeId1]) { indegrees_nodeIdIndex[initNodeId1_index]++; }
				}
			}

			Int knownMinNum = 0;		//已知的排在前面knownMinNum小的indegree
			Int interNodeId_indegree_rank = initSize / 2;		//平均入度的值，也就是所有入度中排在第initSize/2小的入度值，接下来通过快排找到第initSize/2小的入度的initNode的index
			Vec<Index> searchInterNodeIdIndexs(initSize);
			Vec<Index> l_searchInterNodeIdIndexs, r_searchInterNodeIdIndexs;
			for (Int i = 0; i < initSize; ++i) { searchInterNodeIdIndexs[i] = i; }
			l_searchInterNodeIdIndexs.reserve(initSize), r_searchInterNodeIdIndexs.reserve(initSize);

			while (searchInterNodeIdIndexs.size() != 1) {
				l_searchInterNodeIdIndexs.clear(), r_searchInterNodeIdIndexs.clear();
				Index referNodeIdIndex_index = searchInterNodeIdIndexs.size() / 2;
				Index referNodeIdIndex = searchInterNodeIdIndexs[referNodeIdIndex_index];
				Int referNodeIdIndex_indegree = indegrees_nodeIdIndex[referNodeIdIndex];
				
				for (auto searchInterNodeIdIndex = searchInterNodeIdIndexs.begin(); searchInterNodeIdIndex != searchInterNodeIdIndexs.end(); ++searchInterNodeIdIndex) {
					if (searchInterNodeIdIndex == searchInterNodeIdIndexs.begin() + referNodeIdIndex_index) { continue; }

					Int searchInterNodeIdIndex_indegree = indegrees_nodeIdIndex[*searchInterNodeIdIndex];
					if (referNodeIdIndex_indegree >= searchInterNodeIdIndex_indegree) { l_searchInterNodeIdIndexs.emplace_back(*searchInterNodeIdIndex); }
					else { r_searchInterNodeIdIndexs.emplace_back(*searchInterNodeIdIndex); }
				}

				if (knownMinNum + l_searchInterNodeIdIndexs.size() == interNodeId_indegree_rank - 1) {
					searchInterNodeIdIndexs.clear();
					searchInterNodeIdIndexs.emplace_back(referNodeIdIndex);
				}
				else if (knownMinNum + l_searchInterNodeIdIndexs.size() > interNodeId_indegree_rank - 1) {
					searchInterNodeIdIndexs = move(l_searchInterNodeIdIndexs);
				}
				else {
					searchInterNodeIdIndexs = move(r_searchInterNodeIdIndexs);
					knownMinNum += (l_searchInterNodeIdIndexs.size() + 1);
				}
			}

			Index interNodeId_index = searchInterNodeIdIndexs[0];
			NodeId interNodeId = initNodeIds[interNodeId_index];
			Int interNodeId_indegree = indegrees_nodeIdIndex[interNodeId_index];
			Vec<SortStruct<Index, CrossingNum>> initNodeId_indexs_sortByCrossingGap;
			initNodeId_indexs_sortByCrossingGap.reserve(initSize);

			for (Index initNodeId_index = 0; initNodeId_index < initSize; ++initNodeId_index) {
				NodeId initNodeId = initNodeIds[initNodeId_index];
				initNodeId_indexs_sortByCrossingGap.push_back({ initNodeId_index, scc.cm[initNodeId][interNodeId] - scc.cm[interNodeId][initNodeId] });
			}
			sort(initNodeId_indexs_sortByCrossingGap.begin(), initNodeId_indexs_sortByCrossingGap.end());

			Int endsBlockSize = initSize * cfg.initQuickSortEndsRate;
			if (endsBlockSize == 0) { endsBlockSize = 1; }
			NodeIds l_initNodeIds, mid_initNodeIds, r_initNodeIds;
			l_initNodeIds.reserve(endsBlockSize), mid_initNodeIds.reserve(initSize - 2 * endsBlockSize), r_initNodeIds.reserve(endsBlockSize);
			Vec<Index> l_initIndexs, mid_initIndexs, r_initIndexs;
			l_initIndexs.reserve(endsBlockSize), mid_initIndexs.reserve(initSize - 2 * endsBlockSize), r_initIndexs.reserve(endsBlockSize);

			Index initNodeId_sortIndex = 0;
			for (; initNodeId_sortIndex < endsBlockSize; ++initNodeId_sortIndex) {
				l_initNodeIds.emplace_back(initNodeIds[initNodeId_indexs_sortByCrossingGap[initNodeId_sortIndex].id]);
				l_initIndexs.emplace_back(initIndexs[initNodeId_sortIndex]);
			}
			for (; initNodeId_sortIndex < initSize - endsBlockSize; ++initNodeId_sortIndex) {
				mid_initNodeIds.emplace_back(initNodeIds[initNodeId_indexs_sortByCrossingGap[initNodeId_sortIndex].id]);
				mid_initIndexs.emplace_back(initIndexs[initNodeId_sortIndex]);
			}
			for (; initNodeId_sortIndex < initSize; ++initNodeId_sortIndex) {
				r_initNodeIds.emplace_back(initNodeIds[initNodeId_indexs_sortByCrossingGap[initNodeId_sortIndex].id]);
				r_initIndexs.emplace_back(initIndexs[initNodeId_sortIndex]);
			}

			quickSortBlock(l_initNodeIds, l_initIndexs, scc);
			quickSortBlock(r_initNodeIds, r_initIndexs, scc);
			quickSortBlock(mid_initNodeIds, mid_initIndexs, scc);
		}

		void OcmByKq::initQuickSort(SCC& scc) {
			NodeIds initNodeIds(scc.nodeNum);
			Vec<Index> initIndexs(scc.nodeNum);
			for (Int i = 0; i < scc.nodeNum; ++i) {
				initNodeIds[i] = i;
				initIndexs[i] = i;
			}

			quickSortBlock(initNodeIds, initIndexs, scc);

			calculateCN(scc.cm, scc.bestSol);
			return;
		}

		void OcmByKq::mergeSortBlock(NodeIds& sortingNodeIds, Vec<Vec<CrossingNum>>& cm) {
			if (sortingNodeIds.size() == 1) {
				return;
			}

			Int sortingNodeIds_size = sortingNodeIds.size();
			Int l_sortingNodeIds_size = sortingNodeIds_size / 2, r_sortingNodeIds_size = ceil((Real)sortingNodeIds_size / 2);

			NodeIds l_sortingNodeIds(l_sortingNodeIds_size), r_sortingNodeIds(r_sortingNodeIds_size);
			Index pIndex = 0;
			for (Index l_sortingNodeId_index = 0; l_sortingNodeId_index < l_sortingNodeIds_size; ++pIndex, ++l_sortingNodeId_index) { l_sortingNodeIds[l_sortingNodeId_index] = move(sortingNodeIds[pIndex]); }
			for (Index r_sortingNodeId_index = 0; r_sortingNodeId_index < r_sortingNodeIds_size; ++pIndex, ++r_sortingNodeId_index) { r_sortingNodeIds[r_sortingNodeId_index] = move(sortingNodeIds[pIndex]); }
			sortingNodeIds.clear();

			mergeSortBlock(l_sortingNodeIds, cm);
			mergeSortBlock(r_sortingNodeIds, cm);
			//sortingNodeIds.reserve(sortingNodeIds_size);


			if (cfg.mergeSortFlag == Configuration::MergeSortFlag::MergeByClass) {
				Index l_sortingNodeId_index = 0, r_sortingNodeId_index = 0;
				while (l_sortingNodeId_index < l_sortingNodeIds_size && r_sortingNodeId_index < r_sortingNodeIds_size) {
					NodeId l_sortingNodeId = l_sortingNodeIds[l_sortingNodeId_index], r_sortingNodeId = r_sortingNodeIds[r_sortingNodeId_index];
					if (cm[l_sortingNodeId][r_sortingNodeId] <= cm[r_sortingNodeId][l_sortingNodeId]) {
						sortingNodeIds.emplace_back(l_sortingNodeId);
						l_sortingNodeId_index++;
					}
					else {
						sortingNodeIds.emplace_back(r_sortingNodeId);
						r_sortingNodeId_index++;
					}
				}
				while (l_sortingNodeId_index < l_sortingNodeIds_size) { sortingNodeIds.emplace_back(l_sortingNodeIds[l_sortingNodeId_index++]); }
				while (r_sortingNodeId_index < r_sortingNodeIds_size) { sortingNodeIds.emplace_back(r_sortingNodeIds[r_sortingNodeId_index++]); }
			}
			else if (cfg.mergeSortFlag == Configuration::MergeSortFlag::MergeByDP) {
				Vec<Vec<CrossingNum>> crossings_index_pos(r_sortingNodeIds_size, Vec<CrossingNum>(l_sortingNodeIds_size + 1, 0));
				Vec<Vec<CrossingNum>> dp_sumCrossings_index_pos(r_sortingNodeIds_size, Vec<CrossingNum>(l_sortingNodeIds_size + 1, 0));
				Vec<Vec<CrossingNum>> min_dp(r_sortingNodeIds_size, Vec<CrossingNum>(l_sortingNodeIds_size + 1, 0));		//min_dp records[index][pos] min(min_dp[index][0...pos])
				Vec<Vec<Index>> min_dp_poss(r_sortingNodeIds_size, Vec<Index>(l_sortingNodeIds_size + 1, -1));

				for (Index r_sortingNodeId_index = 0; r_sortingNodeId_index < r_sortingNodeIds_size; ++r_sortingNodeId_index) {
					NodeId r_sortingNodeId = r_sortingNodeIds[r_sortingNodeId_index];
					for (Index l_sortingNodeId_index = 0; l_sortingNodeId_index < l_sortingNodeIds_size; ++l_sortingNodeId_index) {
						NodeId l_sortingNodeId = l_sortingNodeIds[l_sortingNodeId_index];
						crossings_index_pos[r_sortingNodeId_index][0] += cm[r_sortingNodeId][l_sortingNodeId];
					}
				}
				for (Index r_sortingNodeId_index = 0; r_sortingNodeId_index < r_sortingNodeIds_size; ++r_sortingNodeId_index) {
					NodeId r_sortingNodeId = r_sortingNodeIds[r_sortingNodeId_index];
					for (Index l_sortingNodeId_index = 0, pos = 1; l_sortingNodeId_index < l_sortingNodeIds_size; ++l_sortingNodeId_index, ++pos) {
						NodeId l_sortingNodeId = l_sortingNodeIds[l_sortingNodeId_index];
						crossings_index_pos[r_sortingNodeId_index][pos] = crossings_index_pos[r_sortingNodeId_index][pos - 1] - cm[r_sortingNodeId][l_sortingNodeId] + cm[l_sortingNodeId][r_sortingNodeId];
					}
				}

				//initialize dp_sumCrossings_index_pos[0]
				dp_sumCrossings_index_pos[0] = move(crossings_index_pos[0]);
				min_dp[0][0] = dp_sumCrossings_index_pos[0][0];
				min_dp_poss[0][0] = 0;
				for (Index pos = 1; pos <= l_sortingNodeIds_size; ++pos) {
					if (dp_sumCrossings_index_pos[0][pos] <= min_dp[0][pos - 1]) {
						min_dp[0][pos] = dp_sumCrossings_index_pos[0][pos];
						min_dp_poss[0][pos] = pos;
					}
					else {
						min_dp[0][pos] = min_dp[0][pos - 1];
						min_dp_poss[0][pos] = min_dp_poss[0][pos - 1];
					}
				}

				//initialize dp_sumCrossings_index_pos[1...r_sortingNodeIs_size][0]
				for (Index r_sortingNodeId_index = 1; r_sortingNodeId_index < r_sortingNodeIds_size; ++r_sortingNodeId_index) {
					dp_sumCrossings_index_pos[r_sortingNodeId_index][0] = crossings_index_pos[r_sortingNodeId_index][0];
					min_dp[r_sortingNodeId_index][0] = crossings_index_pos[r_sortingNodeId_index][0];
					min_dp_poss[r_sortingNodeId_index][0] = 0;
				}

				for (Index r_sortingNodeId_index = 1; r_sortingNodeId_index < r_sortingNodeIds_size; ++r_sortingNodeId_index) {
					for (Index pos = 1; pos <= l_sortingNodeIds_size; ++pos) {
						dp_sumCrossings_index_pos[r_sortingNodeId_index][pos] = min_dp[r_sortingNodeId_index - 1][pos] + crossings_index_pos[r_sortingNodeId_index][pos];

						if (dp_sumCrossings_index_pos[r_sortingNodeId_index][pos] <= min_dp[r_sortingNodeId_index][pos - 1]) {
							min_dp[r_sortingNodeId_index][pos] = dp_sumCrossings_index_pos[r_sortingNodeId_index][pos];
							min_dp_poss[r_sortingNodeId_index][pos] = pos;
						}
						else {
							min_dp[r_sortingNodeId_index][pos] = min_dp[r_sortingNodeId_index][pos - 1];
							min_dp_poss[r_sortingNodeId_index][pos] = pos - 1;
						}
					}
				}

				Vec<Vec<Index>> pos_bucekts(l_sortingNodeIds_size + 1);
				for (Index r_sortingNodeId_index = r_sortingNodeIds_size - 1, pos = l_sortingNodeIds_size; r_sortingNodeId_index >= 0; --r_sortingNodeId_index) {
					Index minPos = min_dp_poss[r_sortingNodeId_index][pos];
					pos_bucekts[minPos].emplace_back(r_sortingNodeId_index);
					pos = minPos;
				}

				for (Index l_sortingNodeId_index = 0; l_sortingNodeId_index < l_sortingNodeIds_size; ++l_sortingNodeId_index) {
					Index pos = l_sortingNodeId_index;
					for (auto r_sortingNodeId_index : pos_bucekts[pos]) { sortingNodeIds.emplace_back(r_sortingNodeIds[r_sortingNodeId_index]); }
					sortingNodeIds.emplace_back(l_sortingNodeIds[l_sortingNodeId_index]);
				}
				for (auto r_sortingNodeId_index : pos_bucekts[l_sortingNodeIds_size]) { sortingNodeIds.emplace_back(r_sortingNodeIds[r_sortingNodeId_index]); }
			}
		}

		void OcmByKq::initMergeSort(Solution& sol, SCC& scc) {
			for (Index i = 0; i < scc.nodeNum; ++i) { scc.bestSol.order[i] = i; }
			//initIndegreeSort(scc);
			mergeSortBlock(sol.order, scc.cm);
			for (Index index = 0; index < scc.nodeNum; ++index) { sol.index[sol.order[index]] = index; }
			calculateCN(scc.cm, sol);
		}

		bool OcmByKq::tarjan(NodeId u, TarjanData& td, Int deepth) {
			if (deepth > cfg.maxTarjanDepth) { return false; }

			td.low[u] = td.dfn[u] = td.dfncnt++, td.sta.push(u), td.inStack[u] = 1;

			for (NodeId v = 0; v < nodeNum; ++v) {
				if (globalCM[u][v] >= globalCM[v][u]) { continue; }

				if (td.dfn[v] == -1) {
					if (!tarjan(v, td, deepth + 1)) { return false; }
					td.low[u] = min(td.low[u], td.low[v]);
				}
				else if (td.inStack[v]) {
					td.low[u] = min(td.low[u], td.dfn[v]);
				}
			}
			if (td.dfn[u] == td.low[u]) {
				while (td.sta.top() != u) {
					NodeId nodeId = td.sta.top();
					td.sta.pop();
					td.scc[td.scccnt].emplace_back(nodeId);
					td.inStack[nodeId] = 0;
				}
				Int nodeId = td.sta.top();
				td.sta.pop();
				td.scc[td.scccnt].emplace_back(nodeId);
				td.inStack[nodeId] = 0;
				td.scccnt++;
			}

			return true;
		}

		void OcmByKq::initSCC(SCC& scc, NodeIds& nodeIds) {
			sort(nodeIds.begin(), nodeIds.end(), less<NodeId>());	// important
			scc.nodeNum = nodeIds.size();
			scc.init(scc.nodeNum, cfg.baseMoveLen, cfg.tabuDividend, cfg.min_tabuMatFixedLen, cfg.min_tabuMatRandLen);
			scc.localNodeId2GlobalNodeId.resize(scc.nodeNum);
			for (Index index = 0; index < scc.nodeNum; ++index) {
				NodeId localNodeId = index;
				scc.localNodeId2GlobalNodeId[localNodeId] = nodeIds[index];
			}

			//initialize scc.cm
			scc.cm.resize(scc.nodeNum);
			for (NodeId localNodeId = 0; localNodeId < scc.nodeNum; ++localNodeId) {
				NodeId globalNodeId = scc.localNodeId2GlobalNodeId[localNodeId];
				scc.cm[localNodeId] = move(globalCM[globalNodeId]);
			}
			for (NodeId localNodeId = 0; localNodeId < scc.nodeNum; ++localNodeId) {
				NodeId globalNodeId = scc.localNodeId2GlobalNodeId[localNodeId];
				for (NodeId nodeIdRow = 0; nodeIdRow < scc.nodeNum; ++nodeIdRow) { scc.cm[nodeIdRow][localNodeId] = scc.cm[nodeIdRow][globalNodeId]; }
			}
		}

		void OcmByKq::initSCCs() {
 			TarjanData td(nodeNum);
			for (NodeId nodeId = 0; nodeId < nodeNum; ++nodeId) {
				if (td.dfn[nodeId] == -1) {
					if (!tarjan(nodeId, td, 0)) {
						td.scccnt = 1;
						for (NodeId pNodeId = 0; pNodeId < nodeNum; ++pNodeId) {
							td.scc[0].emplace_back(pNodeId);
						}
						break;
					}
				}
			}

			// sort scc by toposort
			Vec<Index> sccSortIndexs(td.scccnt);
			Vec<Index> indegree(td.scccnt, 0);
			Vec<Vec<Index>> head(td.scccnt);
			for (Index index1 = 0; index1 < td.scccnt; ++index1) {
				for (Index index2 = index1 + 1; index2 < td.scccnt; ++index2) {

					for (Int i = 0; i < td.scc[index1].size(); ++i) {
						NodeId nodeId1 = td.scc[index1][i];
						for (auto nodeId2 : td.scc[index2]) {
							if (globalCM[nodeId1][nodeId2] < globalCM[nodeId2][nodeId1]) {
								indegree[index2]++;
								head[index1].emplace_back(index2);

								i = td.scc[index1].size();		//这么设置是为了当存在nodeId1和nodeId2有上下关系时，就可以确定两个强连通分量的上下关系。
								break;								//要注意的是，有时候两个强连通分量里的两个点，可能彼此之间没有上下关系，因为谁在上谁在下交点数都一样，同时它们的邻居不一样，比如A点连左边的2、3号点，B点连左边的1、4号点
							}
							else if (globalCM[nodeId1][nodeId2] > globalCM[nodeId2][nodeId1]) {
								indegree[index1]++;
								head[index2].emplace_back(index1);

								i = td.scc[index1].size();
								break;
							}
						}
					}
				}
			}
			stack<Index> sccSta;	//要注意到的是，可能会有入度相等的强连通分量，之前误以为强连通分量之间的入度是严格加1递增的，所以直接用入度赋值
			for (Index sccIndex = 0; sccIndex < td.scccnt; ++sccIndex) {
				if (indegree[sccIndex] == 0) {
					sccSta.push(sccIndex);
				}
			}
			Index pIndex = 0;
			while (!(sccSta.empty())) {
				Index sccIndex = sccSta.top();
				sccSta.pop();
				sccSortIndexs[pIndex++] = sccIndex;

				for (auto childSccIndex : head[sccIndex]) {
					indegree[childSccIndex]--;
					if (indegree[childSccIndex] == 0) {
						sccSta.push(childSccIndex);
					}
				}
			}

			sccs.resize(td.scccnt);
			for (Int i = 0; i < td.scccnt; ++i) {
				Int chosenSCCIndex = sccSortIndexs[i];
				initSCC(sccs[i], td.scc[chosenSCCIndex]);
			}
		}

		void OcmByKq::initCM() {
			globalCM.resize(nodeNum);
			for (auto& globalCM_list : globalCM) { globalCM_list.resize(nodeNum, 0); }
			Vec<Int> dpNodeNum(pInput->nodeNum1);
			
			for (NodeId nodeId1 = 0; nodeId1 < nodeNum - 1; ++nodeId1) {
				if ((pInput->inNodes[nodeId1]).empty()) { continue; }

				fill(dpNodeNum.begin(), dpNodeNum.end(), 0);

				for (auto& nodeId1_inNode : (pInput->inNodes)[nodeId1]) { dpNodeNum[nodeId1_inNode.inNodeId] = nodeId1_inNode.num; }
				for (NodeId pNodeId = 1; pNodeId < (pInput->nodeNum1); ++pNodeId) { dpNodeNum[pNodeId] += dpNodeNum[pNodeId - 1]; }

				for (NodeId nodeId2 = nodeId1 + 1; nodeId2 < nodeNum; ++nodeId2) {
					for (auto nodeId2_inNode : (pInput->inNodes)[nodeId2]) {
						if (nodeId2_inNode.inNodeId != 0) { globalCM[nodeId2][nodeId1] += (dpNodeNum[nodeId2_inNode.inNodeId - 1] * nodeId2_inNode.num); }
						globalCM[nodeId1][nodeId2] += ((dpNodeNum[(pInput->nodeNum1) - 1] - dpNodeNum[nodeId2_inNode.inNodeId]) * nodeId2_inNode.num);
					}
				}
			}

			//Switch to relative CM
			for (NodeId nodeId1 = 0; nodeId1 < nodeNum - 1; ++nodeId1) {
				for (NodeId nodeId2 = nodeId1 + 1; nodeId2 < nodeNum; ++nodeId2) {
					if (globalCM[nodeId1][nodeId2] <= globalCM[nodeId2][nodeId1]) {
						globalCM[nodeId2][nodeId1] -= globalCM[nodeId1][nodeId2];
						globalCM[nodeId1][nodeId2] = 0;
					}
					else {
						globalCM[nodeId1][nodeId2] -= globalCM[nodeId2][nodeId1];
						globalCM[nodeId2][nodeId1] = 0;
					}
				}
			}
		}	//initCM

		void OcmByKq::initSolution() {
			for (auto& scc : sccs) {
				if (scc.nodeNum <= 2) {
					initIndegreeSort(scc);
					scc.isSorted = true;
				}
				/*else if (30 < scc.nodeNum && scc.nodeNum <= 40) {
					//AllNodeBestImp
					scc.ab_conALLNodeIterLen = 2000;
					scc.ab_conBestImpIterLen = 2000;

					for (Index index = 0; index < scc.nodeNum; ++index) {
						scc.bestSol.order[index] = index;
						scc.bestSol.index[index] = index;
					}
					//sortSCCbyExactAlg(scc);
					scc.isSmallSize = true;
					scc.isSorted = false;
					calculateCN(scc.cm, scc.bestSol);
				}*/
				else if (scc.nodeNum <= cfg.SCCSortLen) {
					//DynamicNodesBestImp
					scc.dn_conBestImpIterLen = cfg.dn_smallSCC_conBestImpIterLen;

					switch (cfg.initAlg)
					{
					case Configuration::InitAlgorithm::IndegreeSort:
						initIndegreeSort(scc);
						break;
					case Configuration::InitAlgorithm::Linear:
						for (Index index = 0; index < scc.nodeNum; ++index) {
							scc.bestSol.order[index] = index;
							scc.bestSol.index[index] = index;
						}
						calculateCN(scc.cm, scc.bestSol);
						break;
					case Configuration::InitAlgorithm::CrossingSort :
						initCrossingSort(scc);
						break;
					case Configuration::InitAlgorithm::UpAndDownCrossingSort :
						initUpCrossingAndDownCrossingSort(scc);
						break;
					case Configuration::InitAlgorithm::QuickSort :
						initQuickSort(scc);
						break;
					case Configuration::InitAlgorithm::MergeSort :
						//for (Index i = 0; i < scc.nodeNum; ++i) { scc.bestSol.order[i] = i; }
						initIndegreeSort(scc);
						initMergeSort(scc.bestSol, scc);
						break;
					}
					/*for (Index index = 0; index < scc.nodeNum; ++index) {
						scc.bestSol.order[index] = index;
						scc.bestSol.index[index] = index;
					}*/
					//sortSCCbyExactAlg(scc);
					scc.isSmallSize = true;
					scc.isSorted = false;
					calculateCN(scc.cm, scc.bestSol);
				}
				else {
					//DynamicNodesBestImp
					//scc.dn_conBestImpIterLen = cfg.dn_conBestImpIterLen;
					scc.dn_conBestImpIterLen = scc.nodeNum / cfg.dn_conBestImpIterLenDividend;

					scc.isSorted = false;
					switch (cfg.initAlg)
					{
					case Configuration::InitAlgorithm::IndegreeSort:
						initIndegreeSort(scc);
						break;
					case Configuration::InitAlgorithm::Linear:
						for (Index index = 0; index < scc.nodeNum; ++index) {
							scc.bestSol.order[index] = index;
							scc.bestSol.index[index] = index;
						}
						calculateCN(scc.cm, scc.bestSol);
						break;
					case Configuration::InitAlgorithm::CrossingSort:
						initCrossingSort(scc);
						break;
					case Configuration::InitAlgorithm::UpAndDownCrossingSort :
						initUpCrossingAndDownCrossingSort(scc);
						break;
					case Configuration::InitAlgorithm::QuickSort:
						initQuickSort(scc);
						break;
					case Configuration::InitAlgorithm::MergeSort :
						//for (Index i = 0; i < scc.nodeNum; ++i) { scc.bestSol.order[i] = i; }
						initIndegreeSort(scc);
						initMergeSort(scc.bestSol, scc);
						break;
					}
				}
			}
		}

		void OcmByKq::calculateCN(Vec<Vec<CrossingNum>>& cm, Solution& sol) {
			Int nodeNum = sol.order.size();
			sol.cn = 0;
			for (Index upos = 0; upos < nodeNum; ++upos) {
				NodeId nodeId1 = sol.order[upos];
				for (Index vpos = upos + 1; vpos < nodeNum; ++vpos) {
					NodeId nodeId2 = sol.order[vpos];
					sol.cn += cm[nodeId1][nodeId2];
				}
			}
		}
		#pragma endregion initial

		#pragma region pack

		/*void OcmByKq::pack(SCC& scc) {
			switch (cfg.packFlag)
			{
			case Configuration::PackFlag::PackByRandom :
				packByRandom(scc);
				break;
			case Configuration::PackFlag::PackByRelation :
				packByRelation(scc);
				break;
			}
		}

		void OcmByKq::packByRandom(SCC& scc) {
			scc.packInf.addLayer();
			Int packNodeNum = scc.nodeNum * cfg.packByRandomLowerRate;
			Solution& sol = scc.searchSol;
			Vec<bool> isPackUp(scc.nodeNum, false);		//indexed by nodeId_index
			NodeIds packNodeIds_indexs;
			packNodeIds_indexs.reserve(packNodeNum);
			for (Index packNodeId_index = 1; packNodeId_index <= packNodeNum; ++packNodeId_index) { packNodeIds_indexs.emplace_back(packNodeId_index); }	//因为是往上合并，第0个不能往上合并，所以从1开始索引
			for (Index packNodeId_index = packNodeNum + 1; packNodeId_index < scc.nodeNum; ++packNodeId_index) {
				Int randIndex = rand.pick(packNodeId_index + 1);	//+1是pick是从[0,max)里取随机数
				if (randIndex < packNodeNum) { packNodeIds_indexs[randIndex] = packNodeId_index; }
			}
			for (auto packNodeId_index : packNodeIds_indexs) { isPackUp[packNodeId_index] = true; }

			PackLayerInf& layerInf = scc.packInf.getNowPackLayerInf();
			layerInf.packNodeIdSets[sol.order[0]].emplace_back(sol.order[0]);
			Index prePackNodeId_index = 0;
			NodeId prePackNodeId = sol.order[0];
			for (Index unpackNodeId_index = 1; unpackNodeId_index < scc.nodeNum; ++unpackNodeId_index) {
				if (isPackUp[unpackNodeId_index]) {
					layerInf.packNodeIdSets[prePackNodeId].emplace_back(sol.order[unpackNodeId_index]);
				}
				else {
					sol.order[prePackNodeId_index + 1] = sol.order[unpackNodeId_index];
					sol.index[sol.order[unpackNodeId_index]] = prePackNodeId_index + 1;
					layerInf.packNodeIdSets[sol.order[unpackNodeId_index]].emplace_back(sol.order[unpackNodeId_index]);
					prePackNodeId_index++;
					prePackNodeId = sol.order[unpackNodeId_index];
				}
			}

			scc.nodeNum -= packNodeNum;
			sol.order.resize(scc.nodeNum);
			for (Index index = 0; index < scc.nodeNum; ++index) {
				NodeId packNodeId = sol.order[index];
				NodeIds& nodeIds = layerInf.packNodeIdSets[packNodeId];

				for (auto nodeId1 = nodeIds.begin(); nodeId1 != nodeIds.end(); ++nodeId1) {
					for (auto nodeId2 = nodeId1 + 1; nodeId2 != nodeIds.end(); ++nodeId2) {
						layerInf.sumOfNodeInnerCN += scc.cm[*nodeId1][*nodeId2];
					}
				}
			}
			convertSCC(scc);
			calculateCN(scc.cm, sol);
		}

		void OcmByKq::packByRelation(SCC& scc) {
			scc.packInf.addLayer();
			Int packNodeNum = 0, unsimilarStaNum = cfg.packByRelationUnsimilarStaNum;
			Solution& sol = scc.searchSol;
			PackLayerInf& layerInf = scc.packInf.getNowPackLayerInf();
			Vec<Vec<CrossingNum>>& cm = scc.cm;

			vector<RelationVec> rvs(scc.initNodeSize);
			for (NodeId nodeId = 0; nodeId < scc.initNodeSize; ++nodeId) { rvs[nodeId].init(scc.initNodeSize, nodeId); }
			for (Index index1 = 0; index1 < scc.nodeNum; ++index1) {
				NodeId nodeId1 = sol.order[index1];
				for (Index index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
					NodeId nodeId2 = sol.order[index2];
					if (cm[nodeId1][nodeId2] < cm[nodeId2][nodeId1]) { rvs[nodeId1].setDown(nodeId2); }
					else if (cm[nodeId1][nodeId2] > cm[nodeId2][nodeId1]) { rvs[nodeId2].setDown(nodeId1); }
				}
			}

			for (Index unpackIndex = 0; unpackIndex < scc.nodeNum; ++unpackIndex) {
				bool isSimilar = false;
				NodeId unpackNodeId = sol.order[unpackIndex];
				NodeId similarPackNodeId = -1;

				for (Index packIndex = packNodeNum - 1; packIndex >= 0; --packIndex) {
					NodeId packNodeId = sol.order[packIndex];
					Int unsimilarNum = rvs[packNodeId].countXOR(rvs[unpackNodeId]);

					if (unsimilarNum <= unsimilarStaNum) {
						isSimilar = true;
						similarPackNodeId = packNodeId;
						break;
					}
				}

				if (!isSimilar) {
					NodeId packNodeId = unpackNodeId;
					sol.order[packNodeNum] = unpackNodeId;
					sol.index[unpackNodeId] = packNodeNum;
					layerInf.packNodeIdSets[packNodeId].emplace_back(unpackNodeId);
					packNodeNum++;
				}
				else {
					layerInf.packNodeIdSets[similarPackNodeId].emplace_back(unpackNodeId);
				}
			}

			scc.nodeNum = packNodeNum;
			sol.order.resize(scc.nodeNum);
			for (Index index = 0; index < scc.nodeNum; ++index) {
				NodeId packNodeId = sol.order[index];
				NodeIds& nodeIds = layerInf.packNodeIdSets[packNodeId];

				for (auto nodeId1 = nodeIds.begin(); nodeId1 != nodeIds.end(); ++nodeId1) {
					for (auto nodeId2 = nodeId1 + 1; nodeId2 != nodeIds.end(); ++nodeId2) {
						layerInf.sumOfNodeInnerCN += cm[*nodeId1][*nodeId2];
					}
				}
			}
			convertSCC(scc);
			calculateCN(cm, sol);
		}

		void OcmByKq::unpack(SCC& scc) {
			Solution& sol = scc.searchSol;
			PackLayerInf& layerInf = scc.packInf.getNowPackLayerInf();

			convertBackSCC(scc);

			NodeIds porder = move(sol.order);
			for (auto packNodeId : porder) {
				sol.order.insert(sol.order.end(), layerInf.packNodeIdSets[packNodeId].begin(), layerInf.packNodeIdSets[packNodeId].end());
			}
			Index pindex = 0;
			for (auto unpackNodeId : sol.order) {
				sol.index[unpackNodeId] = pindex++;
			}
			scc.nodeNum = pindex;
			calculateCN(scc.cm, sol);
			scc.packInf.delLayer();
		}

		void OcmByKq::convertSCC(SCC& scc) {
			Solution& sol = scc.searchSol;
			PackLayerInf& layerInf = scc.packInf.getNowPackLayerInf();
			Vec<Vec<CrossingNum>>& cm = scc.cm;

			for (Index packIndex1 = 0; packIndex1 < scc.nodeNum; ++packIndex1) {
				NodeId packNodeId1 = sol.order[packIndex1];
				for (Index packIndex2 = packIndex1 + 1; packIndex2 < scc.nodeNum; ++packIndex2) {
					NodeId packNodeId2 = sol.order[packIndex2];

					CrossingNum pCN1 = 0, pCN2 = 0;
					for (auto unpackNodeId1 : layerInf.packNodeIdSets[packNodeId1]) {
						for (auto unpackNodeId2 : layerInf.packNodeIdSets[packNodeId2]) {
							pCN1 += cm[unpackNodeId1][unpackNodeId2];
							pCN2 += cm[unpackNodeId2][unpackNodeId1];
						}
					}
					cm[packNodeId1][packNodeId2] = pCN1;
					cm[packNodeId2][packNodeId1] = pCN2;
				}
			}
		}

		void OcmByKq::convertBackSCC(SCC& scc) {
			Solution& sol = scc.searchSol;
			PackLayerInf& layerInf = scc.packInf.getNowPackLayerInf();
			Vec<Vec<CrossingNum>>& cm = scc.cm;

			for (Index packIndex1 = 0; packIndex1 < scc.nodeNum; ++packIndex1) {
				NodeId packNodeId1 = sol.order[packIndex1];
				for (Index packIndex2 = packIndex1 + 1; packIndex2 < scc.nodeNum; ++packIndex2) {
					NodeId packNodeId2 = sol.order[packIndex2];

					CrossingNum pCN1 = cm[packNodeId1][packNodeId2], pCN2 = cm[packNodeId2][packNodeId1];
					for (auto unpackNodeId1 : layerInf.packNodeIdSets[packNodeId1]) {
						for (auto unpackNodeId2 : layerInf.packNodeIdSets[packNodeId2]) {
							if (unpackNodeId1 == packNodeId1 && unpackNodeId2 == packNodeId2) { continue; }

							pCN1 -= cm[unpackNodeId1][unpackNodeId2];
							pCN2 -= cm[unpackNodeId2][unpackNodeId1];
						}
					}
					cm[packNodeId1][packNodeId2] = pCN1;
					cm[packNodeId2][packNodeId1] = pCN2;
				}
			}
		}*/
		#pragma endregion pack

		
		/*#pragma region ExactAlg
		void OcmByKq::sortSCCbyGurobi(SCC& scc) {
			MpSolverGurobi mp;
			Vec<Vec<Delta>> w(scc.nodeNum, Vec<Delta>(scc.nodeNum, 0));
			Vec<MpSolverGurobi::Decision> position(scc.nodeNum);
			Vec<Vec<MpSolverGurobi::Decision>> isPrior(scc.nodeNum);

			for (Index index1 = 0; index1 < scc.nodeNum; ++index1) {
				NodeId nodeId1 = scc.bestSol.order[index1];
				for (Index index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
					NodeId nodeId2 = scc.bestSol.order[index2];

					if (scc.cm[nodeId1][nodeId2] < scc.cm[nodeId2][nodeId1]) { w[index1][index2] = scc.cm[nodeId2][nodeId1] - scc.cm[nodeId1][nodeId2]; }
					else if (scc.cm[nodeId1][nodeId2] > scc.cm[nodeId2][nodeId1]) { w[index2][index1] = scc.cm[nodeId1][nodeId2] - scc.cm[nodeId2][nodeId1]; }
				}
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				position[index] = mp.addVar(MpSolverGurobi::VariableType::Integer, 0, scc.nodeNum - 1);
			}

			for (Index index1 = 0; index1 < scc.nodeNum; ++index1) {
				isPrior[index1].resize(scc.nodeNum);
				for (Int index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
					isPrior[index1][index2] = mp.addVar(MpSolverGurobi::VariableType::Bool, 0, 1);
				}
			}

			for (Index index1 = 0; index1 < scc.nodeNum; ++index1) {
				for (Index index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
					MpSolverGurobi::LinearExpr diff = position[index2] - position[index1];
					mp.addConstraint(1 - scc.nodeNum * (1 - isPrior[index1][index2]) <= diff);
					mp.addConstraint(diff <= scc.nodeNum * isPrior[index1][index2] - 1);
				}
			}

			// set objective.
			MpSolverGurobi::LinearExpr objFunc;
			for (Index index1 = 0; index1 < scc.nodeNum; ++index1) {
				for (Index index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
					objFunc += isPrior[index1][index2] * w[index1][index2];
					objFunc += (1 - isPrior[index1][index2]) * w[index2][index1];
				}
			}

			// solve model.
			mp.setOutput(true);
			//mp.setMaxThread(1);
			//mp.setTimeLimit(30 * 60);

			mp.setObjective(objFunc, MpSolverGurobi::OptimaDirection::Maximize);

			// record decision
			Vec<Int> indegree(scc.nodeNum, 0);
			Vec<Vec<Index>> head(scc.nodeNum);
			Vec<Int> order; order.reserve(scc.nodeNum);
			stack<Index> sta;

			if (mp.optimize()) {
				for (Int index1 = 0; index1 < scc.nodeNum; ++index1) {
					for (Int index2 = index1 + 1; index2 < scc.nodeNum; ++index2) {
						if (mp.isTrue(isPrior[index1][index2])) {
							indegree[index2]++;
							head[index1].emplace_back(index2);
						}
						else {
							indegree[index1]++;
							head[index2].emplace_back(index1);
						}
					}
				}
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				if (indegree[index] == 0) { sta.push(index); }
			}

			while (!(sta.empty())) {
				Index indexu = sta.top();
				sta.pop();
				order.emplace_back(indexu);

				for (auto indexv : head[indexu]) {
					indegree[indexv]--;
					if (indegree[indexv] == 0) { sta.push(indexv); }
				}
			}

			Vec<NodeId> pnodeIds(scc.nodeNum);
			for (Index i = 0; i < scc.nodeNum; ++i) {
				pnodeIds[i] = scc.bestSol.order[i];
			}
			for (Int i = 0; i < scc.nodeNum; ++i) {
				Index index = order[i];
				NodeId nodeId = pnodeIds[index];
				Index newIndex = i;

				scc.bestSol.order[newIndex] = nodeId;
				scc.bestSol.index[nodeId] = newIndex;
			}
		}

		void OcmByKq::sortSCCbyExactAlg(SCC& scc) {
			sortSCCbyGurobi(scc);
		}
		#pragma endregion ExactAlg*/

		#pragma region LS
		void OcmByKq::makeMove(Solution& solution, Move& move) {
			if (move.flag == Move::MoveFlag::UP) {
				NodeId moveNodeId = solution.order[move.moveNodeId_index];
				for (Index index = move.moveNodeId_index - 1; index >= move.beMovedNodeId_index; --index) {
					NodeId nodeId = solution.order[index];
					Index new_index = index + 1;
					solution.order[new_index] = nodeId;
					solution.index[nodeId] = new_index;
				}
				solution.order[move.beMovedNodeId_index] = moveNodeId;
				solution.index[moveNodeId] = move.beMovedNodeId_index;
			}
			else if (move.flag == Move::MoveFlag::DOWN) {
				NodeId moveNodeId = solution.order[move.moveNodeId_index];
				for (Index index = move.moveNodeId_index + 1; index <= move.beMovedNodeId_index; ++index) {
					NodeId nodeId = solution.order[index];
					Index new_index = index - 1;
					solution.order[new_index] = nodeId;
					solution.index[nodeId] = new_index;
				}
				solution.order[move.beMovedNodeId_index] = moveNodeId;
				solution.index[moveNodeId] = move.beMovedNodeId_index;
			}

			solution.cn += move.delta;
		}

		void OcmByKq::localSearch_node_once(Int iter, Int l, Int r, Index searchNodeId_index, Move& move, Solution& sol, SCC& scc) {
			NodeId searchNodeId = sol.order[searchNodeId_index];
			Int minDeltaSize = 0;
			Index beMovedNodeId_Index;
			NodeId beMovedNodeId;
			Delta delta, minDelta;

			//UP
			delta = 0;
			for (beMovedNodeId_Index = searchNodeId_index - 1; beMovedNodeId_Index > l; --beMovedNodeId_Index) {
				beMovedNodeId = sol.order[beMovedNodeId_Index];
				delta += (-scc.cm[beMovedNodeId][searchNodeId] + scc.cm[searchNodeId][beMovedNodeId]);

				//TabuMat
				if (beMovedNodeId_Index != 0) {
					if (sol.tabuMat[searchNodeId][beMovedNodeId] >= iter && sol.tabuMat[sol.order[beMovedNodeId_Index - 1]][searchNodeId] >= iter) { continue; }
				}
				else {
					if (sol.tabuMat[searchNodeId][beMovedNodeId] >= iter) { continue; }
				}

				if (minDeltaSize == 0 || minDelta > delta) {
					minDeltaSize = 1;
					minDelta = delta;
					move.flag = Move::MoveFlag::UP;
					move.delta = minDelta;
					move.beMovedNodeId_index = beMovedNodeId_Index;
				}
				else if (minDelta == delta && rand.isPicked(1, ++minDeltaSize)) {
					move.flag = Move::MoveFlag::UP;
					move.beMovedNodeId_index = beMovedNodeId_Index;
				}
			}

			//DOWN
			delta = 0;
			for (beMovedNodeId_Index = searchNodeId_index + 1; beMovedNodeId_Index < r; ++beMovedNodeId_Index) {
				beMovedNodeId = sol.order[beMovedNodeId_Index];
				delta += (-scc.cm[searchNodeId][beMovedNodeId] + scc.cm[beMovedNodeId][searchNodeId]);

				//TabuMat
				if (beMovedNodeId_Index != scc.nodeNum - 1) {
					if (sol.tabuMat[beMovedNodeId][searchNodeId] >= iter && sol.tabuMat[searchNodeId][sol.order[beMovedNodeId_Index + 1]] >= iter) { continue; }
				}
				else {
					if (sol.tabuMat[beMovedNodeId][searchNodeId] >= iter) { continue; }
				}

				if (minDeltaSize == 0 || minDelta > delta) {
					minDeltaSize = 1;
					minDelta = delta;
					move.flag = Move::MoveFlag::DOWN;
					move.delta = minDelta;
					move.beMovedNodeId_index = beMovedNodeId_Index;
				}
				else if (minDelta == delta && rand.isPicked(1, ++minDeltaSize)) {
					move.flag = Move::MoveFlag::DOWN;
					move.beMovedNodeId_index = beMovedNodeId_Index;
				}
			}
		}

		bool OcmByKq::localSearch_BestImp(Int iter, const NodeIds& searchNodeIds, Solution& searchSol, SCC& scc) {
			bool isImproved = false;

			Move bestMove;
			Int bestMove_size = 0;

			for (auto searchNodeId : searchNodeIds) {
				/*if (searchSol.nodeTabuVec[*searchNodeId] > 0) {
					searchSol.nodeTabuVec[*searchNodeId]--;
					continue;
				}*/

				Index searchNodeId_index = searchSol.index[searchNodeId];
				Int l = max(-1, searchNodeId_index - scc.moveLen), r = min(scc.nodeNum, searchNodeId_index + scc.moveLen);
				Move move;
				move.moveNodeId_index = searchNodeId_index;

				localSearch_node_once(iter, l, r, searchNodeId_index, move, searchSol, scc);

				if (move.flag != Move::MoveFlag::NO) {
					if (bestMove_size == 0 || bestMove.delta > move.delta) {
						bestMove = move;
						bestMove_size = 1;
					}
					else if (bestMove.delta == move.delta && rand.isPicked(1, ++bestMove_size)) {
						bestMove = move;
					}
				}
			}

			if (bestMove_size != 0 && bestMove.flag != Move::MoveFlag::NO) {
				if (bestMove.flag == Move::MoveFlag::UP) {
					NodeId moveNodeId = searchSol.order[bestMove.moveNodeId_index], moveNodeId_upNodeId = searchSol.order[bestMove.moveNodeId_index - 1];
					searchSol.tabuMat[moveNodeId_upNodeId][moveNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					if (bestMove.moveNodeId_index != scc.nodeNum - 1) {
						NodeId moveNodeId_downNodeId = searchSol.order[bestMove.moveNodeId_index + 1];
						searchSol.tabuMat[moveNodeId][moveNodeId_downNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					}

					if (bestMove.beMovedNodeId_index != 0) {
						NodeId beMovedNodeId = searchSol.order[bestMove.beMovedNodeId_index], beMovedNodeId_upNodeId = searchSol.order[bestMove.beMovedNodeId_index - 1];
						searchSol.tabuMat[beMovedNodeId_upNodeId][beMovedNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					}
				}
				else {
					NodeId moveNodeId = searchSol.order[bestMove.moveNodeId_index], moveNodeId_downNodeId = searchSol.order[bestMove.moveNodeId_index + 1];
					searchSol.tabuMat[moveNodeId][moveNodeId_downNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					if (bestMove.moveNodeId_index != 0) {
						NodeId moveNodeId_upNodeId = searchSol.order[bestMove.moveNodeId_index - 1];
						searchSol.tabuMat[moveNodeId_upNodeId][moveNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					}

					if (bestMove.beMovedNodeId_index != scc.nodeNum - 1) {
						NodeId beMovedNodeId = searchSol.order[bestMove.beMovedNodeId_index], beMovedNodeId_downNodeId = searchSol.order[bestMove.beMovedNodeId_index + 1];
						searchSol.tabuMat[beMovedNodeId][beMovedNodeId_downNodeId] = iter + scc.tabuMatFixedLen + rand.pick(scc.tabuMatRandLen);
					}
				}

				makeMove(searchSol, bestMove);

				if (scc.bestSol.cn > searchSol.cn) {
					scc.bestSol = scc.searchSol;
					scc.bestSol_size = 1;
					isImproved = true;
				}
				else if (scc.bestSol.cn == searchSol.cn && rand.isPicked(1, ++scc.bestSol_size)) {
					scc.bestSol = scc.searchSol;
					isImproved = true;
				}
			}

			return isImproved;
		}
		#pragma endregion LS

		#pragma region Restart
		void OcmByKq::delByGreedy(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			Int delNodeNum = rand.pick(scc.nodeNum * cfg.minDelNodeNumRate, scc.nodeNum * cfg.maxDelNodeNumRate);
			delNodeIds.reserve(delNodeNum);

			while (delNodeNum--) {
				Vec<CrossingNum> relativeCNs_byNodeId(scc.nodeNum, 0);

				for (Index nodeId1_index = 0; nodeId1_index < sol.order.size(); ++nodeId1_index) {
					NodeId nodeId1 = sol.order[nodeId1_index];
					for (Index nodeId2_index = nodeId1_index + 1; nodeId2_index < sol.order.size(); ++nodeId2_index) {
						NodeId nodeId2 = sol.order[nodeId2_index];

						relativeCNs_byNodeId[nodeId1] += scc.cm[nodeId1][nodeId2];
						relativeCNs_byNodeId[nodeId2] += scc.cm[nodeId1][nodeId2];
					}
				}

				CrossingNum minRelativeCN = -1;
				Index minRelativeCN_index = -1;
				Int minRelativeCN_size = 0;
				for (Index index = 0; index < sol.order.size(); ++index) {
					NodeId nodeId = sol.order[index];

					if (minRelativeCN_size == 0 || minRelativeCN > relativeCNs_byNodeId[nodeId]) {
						minRelativeCN = relativeCNs_byNodeId[nodeId];
						minRelativeCN_index = index;
						minRelativeCN_size = 1;
					}
					else if (minRelativeCN == relativeCNs_byNodeId[nodeId] && rand.isPicked(1, ++minRelativeCN_size)) {
						minRelativeCN_index = index;
					}
				}

				NodeId delNodeId = sol.order[minRelativeCN_index];
				delNodeIds.emplace_back(delNodeId);
				sol.order.erase(sol.order.begin() + minRelativeCN_index);
			}
		}

		void OcmByKq::delByRand(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			Int delNodeNum = rand.pick(scc.nodeNum * cfg.minDelNodeNumRate, scc.nodeNum * cfg.maxDelNodeNumRate);
			delNodeIds.reserve(delNodeNum);
			for (NodeId nodeId = 0; nodeId < delNodeNum; ++nodeId) { delNodeIds.emplace_back(nodeId); }
			for (NodeId nodeId = delNodeNum; nodeId < scc.nodeNum; ++nodeId) {
				Index randIndex = rand.pick(nodeId + 1);	//pick from [0,nodeId + 1)
				if (randIndex < delNodeNum) { delNodeIds[randIndex] = nodeId; }
			}

			Vec<bool> isDel_byNodeId(scc.nodeNum, false);
			for (auto delNodeId : delNodeIds) { isDel_byNodeId[delNodeId] = true; }
			NodeIds pOrder = move(sol.order);
			sol.order.reserve(scc.nodeNum - delNodeNum);
			for (auto nodeId : pOrder) {
				if (!isDel_byNodeId[nodeId]) { sol.order.emplace_back(nodeId); }
			}
		}

		void OcmByKq::delByGreedyRand(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			Int delNodeNum = rand.pick(scc.nodeNum * cfg.minDelNodeNumRate, scc.nodeNum * cfg.maxDelNodeNumRate);
			delNodeIds.reserve(delNodeNum);

			while (delNodeNum--) {
				Vec<CrossingNum> relativeCNs_byNodeId(scc.nodeNum, 0);

				for (Index nodeId1_index = 0; nodeId1_index < sol.order.size(); ++nodeId1_index) {
					NodeId nodeId1 = sol.order[nodeId1_index];
					for (Index nodeId2_index = nodeId1_index + 1; nodeId2_index < sol.order.size(); ++nodeId2_index) {
						NodeId nodeId2 = sol.order[nodeId2_index];

						relativeCNs_byNodeId[nodeId1] += scc.cm[nodeId1][nodeId2];
						relativeCNs_byNodeId[nodeId2] += scc.cm[nodeId1][nodeId2];
					}
				}

				CrossingNum sum_relativeCNs = 0;
				for (auto nodeId : sol.order) { sum_relativeCNs += relativeCNs_byNodeId[nodeId]; }
				if (sum_relativeCNs == 0) {
					NodeId delNodeId_index = rand.pick(sol.order.size());
					NodeId delNodeId = sol.order[delNodeId_index];
					delNodeIds.emplace_back(delNodeId);
					sol.order.erase(sol.order.begin() + delNodeId_index);
				}
				else {
					CrossingNum rand_relativeCNs = rand.pick((CrossingNum)1, sum_relativeCNs + 1);

					NodeId delNodeId = -1;
					Index delNodeId_index = 0;
					CrossingNum pSum_relativeCns = 0;
					for (; delNodeId_index < sol.order.size(); ++delNodeId_index) {
						delNodeId = sol.order[delNodeId_index];
						pSum_relativeCns += relativeCNs_byNodeId[delNodeId];
						if (rand_relativeCNs <= pSum_relativeCns) { break; }
					}

					delNodeIds.emplace_back(delNodeId);
					sol.order.erase(sol.order.begin() + delNodeId_index);
				}
			}
		}

		void OcmByKq::insertByGreedy(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			random_shuffle(delNodeIds.begin(), delNodeIds.end());

			for (auto delNodeId : delNodeIds) {
				CrossingNum relativeCN = 0, minRelativeCN = -1;
				Int minRelativeCN_size = 0;
				Index minRelativeCN_index = -1;

				for (auto orderNodeId : sol.order) { relativeCN += scc.cm[delNodeId][orderNodeId]; }
				minRelativeCN = relativeCN, minRelativeCN_size = 1, minRelativeCN_index = 0;

				for (Index orderNodeId_index = 1; orderNodeId_index <= sol.order.size(); ++orderNodeId_index) {
					NodeId orderNodeId = sol.order[orderNodeId_index - 1];
					relativeCN += (-scc.cm[delNodeId][orderNodeId] + scc.cm[orderNodeId][delNodeId]);

					if (minRelativeCN > relativeCN) {
						minRelativeCN = relativeCN;
						minRelativeCN_index = orderNodeId_index;
						minRelativeCN_size++;
					}
					else if (minRelativeCN == relativeCN && rand.isPicked(1, ++minRelativeCN_size)) {
						minRelativeCN_index = orderNodeId_index;
					}
				}

				sol.order.insert(sol.order.begin() + minRelativeCN_index, delNodeId);
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				sol.index[sol.order[index]] = index;
			}
			calculateCN(scc.cm, sol);
		}

		void OcmByKq::insertByDP(NodeIds& delNodeIds, Solution& sol, SCC& scc) {

		}

		void OcmByKq::insertByRand(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			random_shuffle(delNodeIds.begin(), delNodeIds.end());
			for (auto delNodeId : delNodeIds) {
				Index randIndex = rand.pick(sol.order.size() + 1);
				sol.order.insert(sol.order.begin() + randIndex, delNodeId);
			}
			for (Index index = 0; index < scc.nodeNum; ++index) {
				sol.index[sol.order[index]] = index;
			}
			calculateCN(scc.cm, sol);
		}

		void OcmByKq::insertByGreedyRand(NodeIds& delNodeIds, Solution& sol, SCC& scc) {
			random_shuffle(delNodeIds.begin(), delNodeIds.end());

			for (auto delNodeId : delNodeIds) {
				Vec<CrossingNum> relaticeCNs_byIndex(sol.order.size() + 1);

				for (auto orderNodeId : sol.order) { relaticeCNs_byIndex[0] += scc.cm[delNodeId][orderNodeId]; }

				for (Index orderNodeId_index = 1; orderNodeId_index <= sol.order.size(); ++orderNodeId_index) {
					NodeId orderNodeId = sol.order[orderNodeId_index - 1];
					relaticeCNs_byIndex[orderNodeId_index] += (-scc.cm[delNodeId][orderNodeId] + scc.cm[orderNodeId][delNodeId]);
				}

				CrossingNum sum_relativeCNs = 0;
				for (Index index = 0; index <= sol.order.size(); ++index) { sum_relativeCNs += relaticeCNs_byIndex[index]; }
				if (sum_relativeCNs == 0) {
					Index insertIndex = rand.pick(sol.order.size() + 1);
					sol.order.insert(sol.order.begin() + insertIndex, delNodeId);
				}
				else {
					CrossingNum rand_relativeCNs = rand.pick((CrossingNum)1, sum_relativeCNs + 1);

					Index insertIndex = 0;
					CrossingNum pSum_relativeCNs = 0;
					for (; insertIndex <= sol.order.size(); ++insertIndex) {
						pSum_relativeCNs += relaticeCNs_byIndex[insertIndex];
						if (rand_relativeCNs <= pSum_relativeCNs) { break; }
					}

					sol.order.insert(sol.order.begin() + insertIndex, delNodeId);
				}
			}

			for (Index index = 0; index < scc.nodeNum; ++index) {
				sol.index[sol.order[index]] = index;
			}
			calculateCN(scc.cm, sol);
		}

		void OcmByKq::restartByGreedyAndRand(Solution& sol, SCC& scc) {
			NodeIds delNodeIds;

			switch (cfg.restartByDelFalg)
			{
			case Configuration::RestartByDelFlag::RestartByGreedyDel :
				delByGreedy(delNodeIds, scc.searchSol, scc);
				break;
			case Configuration::RestartByDelFlag::RestartByRandDel :
				delByRand(delNodeIds, scc.searchSol, scc);
				break;
			case Configuration::RestartByDelFlag::RestartByGreedyRandDel :
				delByGreedyRand(delNodeIds, scc.searchSol, scc);
				break;
			}

			switch (cfg.restartByInsertFlag)
			{
			case Configuration::RestartByInsertFlag::RestartByGreedyInsert :
				insertByGreedy(delNodeIds, scc.searchSol, scc);
				break;
			case Configuration::RestartByInsertFlag::RestartByRandInsert :
				insertByRand(delNodeIds, scc.searchSol, scc);
				break;
			case Configuration::RestartByInsertFlag::RestartByGreedyRandInsert :
				insertByGreedyRand(delNodeIds, scc.searchSol, scc);
				break;
			}
		}
		#pragma endregion Restart

		#pragma region HEA
		void OcmByKq::CrossOverByRandom(Solution& sol1, Solution& sol2, Solution& childSol, Vec<Vec<CrossingNum>>& cm) {
			//make sure sol1.order.size = sol2.order.size

			Int nodeNum = sol1.order.size();
			Vec<Degree> indegree(nodeNum, 0);
			Vec<NodeIds> head(nodeNum);
			for (Index index1 = 0; index1 < nodeNum; ++index1) {
				NodeId nodeId1 = sol1.order[index1];
				for (Index index2 = index1 + 1; index2 < nodeNum; ++index2) {
					NodeId nodeId2 = sol2.order[index2];

					if (sol2.index[nodeId1] < sol2.index[nodeId2]) {
						indegree[nodeId2]++;
						head[nodeId1].emplace_back(nodeId2);
					}
				}
			}

			Vec<Index> pindex(nodeNum);
			Vec<NodeIds> porder(nodeNum);
			queue<NodeId> que;
			for (NodeId nodeId = 0; nodeId < nodeNum; ++nodeId) {
				if (indegree[nodeId] == 0) {
					pindex[nodeId] = 0;
					porder[0].emplace_back(nodeId);
					que.push(nodeId);
				}
			}
			while (!(que.empty())) {
				NodeId pnodeId = que.front();
				que.pop();

				for (auto nextNodeId : head[pnodeId]) {
					indegree[nextNodeId]--;
					if (indegree[nextNodeId] == 0) {
						pindex[nextNodeId] = pindex[pnodeId] + 1;
						porder[pindex[nextNodeId]].emplace_back(nextNodeId);
						que.push(nextNodeId);
					}
				}
			}

			for (Index index = 0; index < nodeNum && !(porder[index].empty()); ++index) {
				srand(rand.pick(1000));
				random_shuffle(porder[index].begin(), porder[index].end());
			}

			Index solIndex = 0;
			for (Index index = 0; index < nodeNum && !(porder[index].empty()); ++index) {
				for (auto nodeId : porder[index]) {
					childSol.order[solIndex] = nodeId;
					childSol.index[nodeId] = solIndex;
					solIndex++;
				}
			}
			calculateCN(cm, childSol);
		}

		void OcmByKq::CrossOverByCrossingSort(Solution& sol1, Solution& sol2, Solution& childSol, Vec<Vec<CrossingNum>>& cm) {
			//make sure sol1.order.size = sol2.order.size

			Int nodeNum = sol1.order.size();
			Vec<Degree> indegree(nodeNum, 0);
			Vec<NodeIds> head(nodeNum);
			Vec<CrossingNum> crossings(nodeNum, 0);
			for (Index index1 = 0; index1 < nodeNum; ++index1) {
				NodeId nodeId1 = sol1.order[index1];
				for (Index index2 = index1 + 1; index2 < nodeNum; ++index2) {
					NodeId nodeId2 = sol2.order[index2];

					if (sol2.index[nodeId1] < sol2.index[nodeId2]) {
						indegree[nodeId2]++;
						head[nodeId1].emplace_back(nodeId2);
					}
				}
			}

			ConsecutiveIdSet indegree0_nodeIds(nodeNum);
			for (NodeId nodeId = 0; nodeId < nodeNum; ++nodeId) {
				if (indegree[nodeId] == 0) {
					indegree0_nodeIds.insert(nodeId);
				}
			}
			for (Index indegree0_nodeId1_index = 0; indegree0_nodeId1_index < indegree0_nodeIds.size(); ++indegree0_nodeId1_index) {
				NodeId indegree0_nodeId1 = indegree0_nodeIds.itemAt(indegree0_nodeId1_index);
				for (Index indegree0_nodeId2_index = indegree0_nodeId1_index + 1; indegree0_nodeId2_index < indegree0_nodeIds.size(); ++indegree0_nodeId2_index) {
					NodeId indegree0_nodeId2 = indegree0_nodeIds.itemAt(indegree0_nodeId2_index);

					crossings[indegree0_nodeId1] += cm[indegree0_nodeId1][indegree0_nodeId2];
					crossings[indegree0_nodeId2] += cm[indegree0_nodeId2][indegree0_nodeId1];
				}
			}

			Index sortIndex = 0;
			while (!(indegree0_nodeIds.empty())) {
				//find the nodeId with the min crossing value
				CrossingNum minCro;
				NodeId minCro_nodeId = -1;
				Int minCro_size = 0;
				for (Index indegree0_nodeId_index = 0; indegree0_nodeId_index < indegree0_nodeIds.size(); ++indegree0_nodeId_index) {
					NodeId indegree0_nodeId = indegree0_nodeIds.itemAt(indegree0_nodeId_index);

					if (minCro_size == 0 || minCro > crossings[indegree0_nodeId]) {
						minCro = crossings[indegree0_nodeId];
						minCro_nodeId = indegree0_nodeId;
						minCro_size = 1;
					}
					else if (minCro == crossings[indegree0_nodeId] && rand.isPicked(1, ++minCro_size)) {
						minCro_nodeId = indegree0_nodeId;
					}
				}

				childSol.order[sortIndex] = minCro_nodeId;
				childSol.index[minCro_nodeId] = sortIndex;
				sortIndex++;
				indegree0_nodeIds.eraseItem(minCro_nodeId);
				for (Index indeg0_nodeId_index = 0; indeg0_nodeId_index < indegree0_nodeIds.size(); ++indeg0_nodeId_index) {
					NodeId indeg0_nodeId = indegree0_nodeIds.itemAt(indeg0_nodeId_index);
					crossings[indeg0_nodeId] -= cm[indeg0_nodeId][minCro_nodeId];
				}

				for (auto nextNodeId : head[minCro_nodeId]) {
					indegree[nextNodeId]--;
					if (indegree[nextNodeId] == 0) {
						for (Index indeg0_nodeId_index = 0; indeg0_nodeId_index < indegree0_nodeIds.size(); ++indeg0_nodeId_index) {
							NodeId indeg0_nodeId = indegree0_nodeIds.itemAt(indeg0_nodeId_index);
							crossings[indeg0_nodeId] += cm[indeg0_nodeId][nextNodeId];
							crossings[nextNodeId] += cm[nextNodeId][indeg0_nodeId];
						}

						indegree0_nodeIds.insert(nextNodeId);
					}
				}
			}
			calculateCN(cm, childSol);
		}
		#pragma endregion HEA
		OcmByKq::Performance OcmByKq::solve_DynamicNodes_BestImp_LocalSearch(Output& output) {
			OcmByKq::Performance perf;

			Int piter = 0;
			Vec<Int> stageIters(sccs.size(), 0);
			Vec<Int> conBestImpIters(sccs.size(), 0);
			Vec<Int> dynamicNodesSizes(sccs.size(), 1);
			Vec<Int> packDeepths(sccs.size(), 0);
			for (Index index = 0; index < sccs.size(); ++index) {
				if (sccs[index].isSorted) { continue; }
				sccs[index].searchSol = sccs[index].bestSol;
			}

			bool whileFlag = true;
			while (whileFlag) {
				whileFlag = false;
				for (Index sccIndex = 0; !isBreakOut() && sccIndex < sccs.size(); ++sccIndex) {
					if (sccs[sccIndex].isSorted) { continue; }
					whileFlag = true;

					Int& stageIter = stageIters[sccIndex];
					Int& conBestImpIter = conBestImpIters[sccIndex];
					Int& dynamicNodesSize = dynamicNodesSizes[sccIndex];
					Int& packDeepth = packDeepths[sccIndex];
					SCC& scc = sccs[sccIndex];
					Solution& searchSol = scc.searchSol;

					Int searchDynamicNodeIds_vec_size = -1;
					Vec<NodeIds> searchDynamicNodeIds_vec;
					//Vec<Int> searchOrder;
					if (scc.dynamicNodesBestImpLsFlag == SCC::DynamicNodesBestImpLSFlag::dnBestImp) {
						searchDynamicNodeIds_vec_size = ceil((Real)(scc.nodeNum) / dynamicNodesSize);
						searchDynamicNodeIds_vec.resize(searchDynamicNodeIds_vec_size);
						for (auto& searchDynamicNodeIds : searchDynamicNodeIds_vec) { searchDynamicNodeIds.reserve(dynamicNodesSize); }

						Index pIndex = 0;
						for (Index searchDynamicNodeIds_index = 0; searchDynamicNodeIds_index < searchDynamicNodeIds_vec_size; ++searchDynamicNodeIds_index) {
							Int pSize = dynamicNodesSize;
							while (pSize-- && pIndex < scc.nodeNum) {
								NodeId searchNodeId = scc.searchSol.order[pIndex++];
								searchDynamicNodeIds_vec[searchDynamicNodeIds_index].emplace_back(searchNodeId);
							}
						}

						//searchOrder.reserve(searchDynamicNodeIds_vec_size);
						//for (Int i = 0; i < searchDynamicNodeIds_vec_size; ++i) { searchOrder.emplace_back(i); }
					}

					bool isImproved = false;
					switch (scc.dynamicNodesBestImpLsFlag)
					{
					case SCC::DynamicNodesBestImpLSFlag::dnBestImp :
						/*if (dynamicNodesSize > cfg.dn_startNodeTabuBound) {
							minNodeTabuLen = dynamicNodesSize / 100;
							for (auto nodeTabuBound : cfg.dn_nodeTabuBounds) {
								if (dynamicNodesSize <= nodeTabuBound) { minNodeTabuLen = nodeTabuBound / 100; }
							}
						}*/

						for (auto& searchDynamicNodeIds : searchDynamicNodeIds_vec) {
							if (localSearch_BestImp(stageIter, searchDynamicNodeIds, searchSol, scc)) {
								isImproved = true;
							}
						}

						//disorder
						/*random_shuffle(searchOrder.begin(), searchOrder.end());
						for (Index searchDynamicNodeIds_index : searchOrder) {
							if (localSearch_BestImp(stageIter, 1, searchDynamicNodeIds_vec[searchDynamicNodeIds_index], searchSol, scc)) {
								isImproved = true;
							}
						}*/

						if (isImproved) {
							scc.moveLen -= cfg.dn_MoveLenLowerGap;
							scc.moveLen = max(min(scc.nodeNum, cfg.baseMoveLen * dynamicNodesSize), scc.moveLen);
							conBestImpIter = 0;
						}
						else {
							if (conBestImpIter == scc.dn_conBestImpIterLen) {
								if (dynamicNodesSize == scc.nodeNum) {
									if (scc.isSmallSize) { scc.isSorted = true; }
									else { scc.dynamicNodesBestImpLsFlag = SCC::DynamicNodesBestImpLSFlag::dnRestart; }
								}
								else {
									dynamicNodesSize = min(scc.nodeNum, dynamicNodesSize * cfg.dn_dynamicNodeIdsSizeUpperRate);
									scc.moveLen = min(scc.nodeNum, dynamicNodesSize * cfg.baseMoveLen);
									//scc.searchSol.clearNodeTabuVec();
									scc.searchSol.clearTabuMat();
									stageIter = 0;
									conBestImpIter = 0;
								}
							}
							else {
								/*if (scc.moveLen == scc.nodeNum) { conBestImpIter++; }
								else {
									scc.moveLen *= cfg.dn_MoveLenUpperGap;
									scc.moveLen = min(scc.nodeNum, scc.moveLen);
								}*/
								conBestImpIter++;
								scc.moveLen *= cfg.dn_MoveLenUpperGap;
								scc.moveLen = min(scc.nodeNum, scc.moveLen);
							}
						}
						stageIter++;
						break;
					case SCC::DynamicNodesBestImpLSFlag::dnRestart :
						switch (cfg.restartAlgFlag)
						{
						case Configuration::RestartAlgFlag::RestartByRandomCrossOver :
							CrossOverByRandom(scc.bestSol, scc.searchSol, scc.searchSol, scc.cm);
							break;
						case Configuration::RestartAlgFlag::RestartByGreedyCrossOver :
							CrossOverByCrossingSort(scc.bestSol, scc.searchSol, scc.searchSol, scc.cm);
							break;
						case Configuration::RestartAlgFlag::RestartByMergeSort :
							//cfg.mergeSortFlag = Configuration::MergeSortFlag::MergeByClass;
							initMergeSort(scc.searchSol, scc);
							break;
						case Configuration::RestartAlgFlag::RestartByGreedyAndRand :
							restartByGreedyAndRand(scc.searchSol, scc);
							break;
						}
						scc.dynamicNodesBestImpLsFlag = SCC::DynamicNodesBestImpLSFlag::dnBestImp;
						dynamicNodesSize = 1;
						scc.moveLen = cfg.baseMoveLen;
						scc.searchSol.clearTabuMat();
						stageIter = 0;
						conBestImpIter = 0;
						break;
					}
				}
				piter++;

				signal(SIGTERM, sigtermHandler);
			}

			//save output
			saveOutput(sccs, output);

			perf.stopState += ("iter = " + to_string(piter));
			perf.msConvergenceCpu = timer.elapsedMilliseconds();
			return perf;
		}
	}	//namespace OCM
}	//namespace GOAL