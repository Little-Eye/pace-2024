////////////////////////////////
/// usage : 1.	ocm algorithm.
///			2.	pace2023
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_OCM_BY_KQ_H
#define CN_HUST_GOAL_OCM_BY_KQ_H


#include <algorithm>
#include <queue>
#include <math.h>

#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Common/Arr.h"
#include "GOAL/Common/Log.h"
#include "GOAL/Common/Timer.h"
#include "GOAL/Common/Random.h"
#include "GOAL/Common/ConsecutiveIdMap.h"
#include "GOAL/Common/ConsecutiveIdSet.h"
#include "GOAL/Optimization/Interface/ConfigurationBase.h"
#include "GOAL/Optimization/Interface/PerformanceBase.h"
#include "GOAL/Optimization/OCM/Problem.h"
#include "GOAL/Optimization/OCM/Visualization.h"
//#include "GOAL/ThirdParty/Gurobi.h"


namespace goal {
	namespace OCM {
		class OcmByKq {
		public:
			struct Environment : public EnvBase {};

			struct Configuration : public CfgBase {
				enum Algorithm {
					DynamicNodesBestImpLocalSearch
				};

				enum InitAlgorithm {
					IndegreeSort,
					Linear,
					CrossingSort,
					UpAndDownCrossingSort,
					QuickSort,
					MergeSort
				};

				enum NodeTabuRateMinusFlag {
					nodeTabuRateMinusByGap,
					nodeTabuRateMinusBy0
				};

				enum RestartAlgFlag {
					RestartByRandomCrossOver,
					RestartByGreedyCrossOver,
					RestartByMergeSort,
					RestartByGreedyAndRand
				};

				enum MergeSortFlag {
					MergeByClass,
					MergeByDP
				};

				enum RestartByDelFlag {
					RestartByGreedyDel,
					RestartByRandDel,
					RestartByGreedyRandDel
				};

				enum RestartByInsertFlag {
					RestartByGreedyInsert,
					RestartByRandInsert,
					RestartByGreedyRandInsert
				};

				virtual Str briefStr() const override;

				Algorithm alg = Algorithm::DynamicNodesBestImpLocalSearch;
				InitAlgorithm initAlg = InitAlgorithm::IndegreeSort;
				MergeSortFlag mergeSortFlag = MergeSortFlag::MergeByDP;
				RestartAlgFlag restartAlgFlag = RestartAlgFlag::RestartByGreedyAndRand;
				RestartByDelFlag restartByDelFalg = RestartByDelFlag::RestartByGreedyRandDel;
				RestartByInsertFlag restartByInsertFlag = RestartByInsertFlag::RestartByGreedyRandInsert;

				Int maxTarjanDepth = 5000;
				Int SCCSortLen = 50;
				Int baseMoveLen = 10;

				//Tabu
				Int tabuDividend = 100;
				Int min_tabuMatFixedLen = 3;
				Int min_tabuMatRandLen = 3;

				//DynamicNodesBestImp
				Int dn_smallSCC_conBestImpIterLen = 2000;
				//Int dn_conBestImpIterLen = 100;
				Int dn_conBestImpIterLenDividend = 7;
				Int dn_MoveLenUpperGap = 2;
				Int dn_MoveLenLowerGap = 1;
				Int dn_dynamicNodeIdsSizeUpperRate = 2;

				//init
				Real initQuickSortEndsRate = 0.2;		//快排求初解时每次分三段，此参数控制前段和后段节点所占总的比例，必须小于0.5

				//Restart
				Real minDelNodeNumRate = 0.3;
				Real maxDelNodeNumRate = 0.5;
			};

			struct Performance : public PerfWhiteBoxBase {};


			static constexpr ID HoleIdOfBoundary = -1;
			static constexpr ID NoId = -1;

			bool init(const Input& input, const Environment& env,
				const Configuration& cfg,
				const Record<CrossingNum>& rec = Record<CrossingNum>());

			Performance solve(Output& output);

			#pragma region Algorithm
			bool isBreakOut();
			void saveOutput(Vec<SCC>& sccs, Output& output);
			void calculateCN(Vec<Vec<CrossingNum>>& cm, Solution& sol);
			bool tarjan(NodeId u, TarjanData& td, Int deepth);
			//void sortSCCbyGurobi(SCC& scc);
			//void sortSCCbyExactAlg(SCC& scc);

			//LS
			void makeMove(Solution& solution, Move& move);
			//void localSearch_node_once(Int iter, Int l, Int r, SCC& scc, Solution& sol, Move& move, Index indexu);
			void localSearch_node_once(Int iter, Int l, Int r, Index searchNodeId_index, Move& move, Solution& sol, SCC& scc);
			//bool localSearch_BestImp(Int iter, const NodeIds& searchNodeIds, Solution& searchSol, SCC& scc);
			bool localSearch_BestImp(Int iter, const NodeIds& searchNodeIds, Solution& searchSol, SCC& scc);

			//HEA
			//void updatePool(Population& pop);
			void CrossOverByRandom(Solution& sol1, Solution& sol2, Solution& childSol, Vec<Vec<CrossingNum>>& cm);
			void CrossOverByCrossingSort(Solution& sol1, Solution& sol2, Solution& childSol, Vec<Vec<CrossingNum>>& cm);
			//Solution& getChildByCrossover(SCC& scc, Population& pop);

			//Performance solvePartTabuBasedLoaclSearch(Output& output);
			Performance solve_DynamicNodes_BestImp_LocalSearch(Output& output);
			//Performance solveCountLowerBound(Output& output);
			//Performance solveHybridEvolutionAlgrithm(Output& output);
			//Performance solveGreedy(Output& output);
			//Performance solveLocalSearch(Output& output);
			#pragma endregion Algorithm

			#pragma region DataStructure
			#pragma endregion DataStructure

			#pragma region Utility
			const Input& input() const { return *pInput; }
			#pragma endregion Utility

			#pragma region Preprocess
			void initCM();
			void initSCCs();
			void initSolution();
			void initSCC(SCC& scc, NodeIds& nodeIds);
			void initIndegreeSort(SCC& scc);
			void initCrossingSort(SCC& scc);
			void initUpCrossingAndDownCrossingSort(SCC& scc);
			void initQuickSort(SCC& scc);
			void quickSortBlock(NodeIds& initNodeIds, Vec<Index>& initIndexs, SCC& scc);
			void initMergeSort(Solution& sol, SCC& scc);
			void mergeSortBlock(NodeIds& sortingNodeIds, Vec<Vec<CrossingNum>>& cm);
			#pragma endregion Preprocess

			#pragma region PackProcess
			/*void pack(SCC& scc);
			void unpack(SCC& scc);
			void packByRandom(SCC& scc);
			void packByRelation(SCC& scc);
			void convertSCC(SCC& scc);
			void convertBackSCC(SCC& scc);*/
			#pragma endregion PackProcess

			#pragma region Restart
			void delByGreedy(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void delByRand(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void delByGreedyRand(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void insertByGreedy(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void insertByDP(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void insertByRand(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void insertByGreedyRand(NodeIds& delNodeIds, Solution& sol, SCC& scc);
			void restartByGreedyAndRand(Solution& sol, SCC& scc);
			#pragma region Restart

			#pragma region Tabu
			#pragma endregion Tabu
		protected:
			const Input* pInput;

			Environment env;
			Configuration cfg;
			Random rand;
			Timer timer;

			#pragma warning(push)
			#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
			Int nodeNum;
			Vec<Vec<CrossingNum>> globalCM;		//crossing matrix
			Vec<SCC> sccs;

			#pragma warning(pop)

			#pragma region Solve
			#pragma endregion Solve
		};

	}
}


#endif // CN_HUST_GOAL_OCM_BY_KQ_H
