////////////////////////////////
/// usage : 1.	strongly connected component for directed graph.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_GRAPH_STRONGLY_CONNECTED_COMPONENT_H
#define CN_HUST_GOAL_COMMON_GRAPH_STRONGLY_CONNECTED_COMPONENT_H


#include "GOAL/Typedef.h"
#include "GOAL/Common/Arr.h"
#include "GOAL/Common/Math.h"
#include "GOAL/Common/Container.h"
#include "GOAL/Graph/Graph.h"
#include "GOAL/Graph/DepthFirstSearch.h"


namespace goal {
namespace graph {

template<typename AdjNode = ID>
class StronglyConnectedComponentR { // recursive implementation.
public:
    using AdjNodes = graph::AdjNodes<AdjNode>;


    static constexpr ID InvalidId = -1;


    void solve(const AdjList<AdjNode>& adjacencyList) {
        ID nodeNum = sCast<ID>(adjacencyList.size());
        adjList = &adjacencyList;
        ages.clear();
        ages.resize(nodeNum, InvalidId);
        oldestReachableAncestors.clear();
        oldestReachableAncestors.resize(nodeNum, InvalidId);
        inStack.clear();
        inStack.resize(nodeNum, false);
        nodeStack.clear();
        nodeStack.reserve(nodeNum);
        scc.resize(nodeNum);

        age = 0;
        for (ID i = 0; i < nodeNum; ++i) { if (ages[i] == InvalidId) { expand(i); } }
    }

    const AdjNodes& adjNodes(ID n) const { return (*adjList)[n]; }

    void expand(ID r) {
        ages[r] = oldestReachableAncestors[r] = age++;
        nodeStack.push_back(r);
        inStack[r] = true;

        for (auto n = adjNodes(r).begin(); n != adjNodes(r).end(); ++n) { // for each neighbor.
            if (ages[*n] == InvalidId) { // if not visited, check if the subtree rooted with `*n` connects to ancestors of `r`.
                expand(*n);
                math::updateMin(oldestReachableAncestors[r], oldestReachableAncestors[*n]);
			} else if (inStack[*n]) { // `*n` is still in stack, i.e., it is a back edge, not cross edge.
                math::updateMin(oldestReachableAncestors[r], ages[*n]);
            }
        }

        if (oldestReachableAncestors[r] == ages[r]) { // head node found, pop the stack.
            ID sscId = sCast<ID>(sccs.size());
            Vec<ID>& component(newBack(sccs));
            ID t;
            do {
                t = nodeStack.back();
                component.push_back(t);
                scc[t] = sscId;
                inStack[t] = false;
                nodeStack.pop_back();
            } while (t != r);
        }
    }

    // input.
    const AdjList<AdjNode>* adjList;

    // temp data.
    ID age; // smaller is older.
    Vec<ID> ages; // `ages[n]` is the sequence of visit (topo order) of node `n`.
    Vec<ID> oldestReachableAncestors; // `oldestReachableAncestors[n]` is the earliest visited node reachable from subtree rooted with node `n`.
    Vec<bool> inStack;
    Vec<ID> nodeStack; // all the connected ancestors (could be part of SCC).

    // output.
    Vec<Vec<ID>> sccs; // `sccs[k]` is the nodes in the `k`th SCC.
    Vec<ID> scc; // `scc[n]` is the SCC of node `n`.
};

template<typename AdjNode = ID>
class StronglyConnectedComponentNR { // non-recursive implementation.
public:
    static constexpr ID InvalidId = -1;


    void solve(const AdjList<AdjNode>& adjList) {
        ID nodeNum = sCast<ID>(adjList.size());
        ages.clear();
        ages.resize(nodeNum, InvalidId);
        oldestReachableAncestors.clear();
        oldestReachableAncestors.resize(nodeNum, InvalidId);
        inStack.clear();
        inStack.resize(nodeNum, false);
        nodeStack.clear();
        nodeStack.reserve(nodeNum);
        scc.resize(nodeNum);

		age = 0;
        Vec<ID> neighborIters;
		for (ID i = 0; i < nodeNum; ++i) {
			if (ages[i] != InvalidId) { continue; }

            auto visitNodeOnPush = [&](ID r, ID) ->bool {
                ages[r] = oldestReachableAncestors[r] = age++;
                nodeStack.push_back(r);
                inStack[r] = true;
                return false;
            };

            auto visitNodeBeforePop = [&](ID r) ->bool {
                if (oldestReachableAncestors[r] == ages[r]) { // head node found, pop the stack.
                    ID sscId = sCast<ID>(sccs.size());
                    Vec<ID>& component(newBack(sccs));
                    ID t;
                    do {
                        t = nodeStack.back();
                        component.push_back(t);
                        scc[t] = sscId;
                        inStack[t] = false;
                        nodeStack.pop_back();
                    } while (t != r);
                }
                return false;
            };
            auto visitNodeOnPop = [&](ID n, ID r) ->bool {
                visitNodeBeforePop(n);
                math::updateMin(oldestReachableAncestors[r], oldestReachableAncestors[n]);
                return false;
            };

            neighborIters.clear();
            neighborIters.resize(nodeNum, 0);
            auto nextNeighbor = [&](ID from) ->ID {
                if (neighborIters[from] >= adjList[from].size()) { return InvalidId; }
                return adjList[from][neighborIters[from]++];
            };

            auto isVisitedNode = [&](ID n, ID r) ->bool {
				if (ages[n] == InvalidId) { return false; }
				if (inStack[n]) { math::updateMin(oldestReachableAncestors[r], ages[n]); }
                return true;
            };

            visitNodeOnPush(i, InvalidId);
			dfs::depthFirstSearch<ID>(nodeNum, i, visitNodeOnPush, visitNodeOnPop, nextNeighbor, isVisitedNode);
            visitNodeBeforePop(i);
		}
	}

    // temp data.
    ID age; // smaller is older.
    Vec<ID> ages; // `ages[n]` is the sequence of visit (topo order) of node `n`.
    Vec<ID> oldestReachableAncestors; // `oldestReachableAncestors[n]` is the earliest visited node reachable from subtree rooted with node `n`.
    Vec<bool> inStack;
    Vec<ID> nodeStack; // all the connected ancestors (could be part of SCC).

    // output.
    Vec<Vec<ID>> sccs; // `sccs[k]` is the nodes in the `k`th SCC.
    Vec<ID> scc; // `scc[n]` is the SCC of node `n`.
};

template<typename AdjNode = ID>
using StronglyConnectedComponent = StronglyConnectedComponentNR<AdjNode>;

}
}


#endif // CN_HUST_GOAL_COMMON_GRAPH_STRONGLY_CONNECTED_COMPONENT_H

// sample code.
//ID n = 9;
//graph::AdjList<ID> adjList(n);
//adjList[0].resize(1);
//adjList[0][0] = 1;
//
//adjList[1].resize(2);
//adjList[1][0] = 2;
//adjList[1][1] = 5;
//
//adjList[2].resize(2);
//adjList[2][0] = 3;
//adjList[2][1] = 4;
//
//adjList[3].resize(1);
//adjList[3][0] = 1;
//
//adjList[5].resize(2);
//adjList[5][0] = 4;
//adjList[5][1] = 6;
//
//adjList[6].resize(1);
//adjList[6][0] = 7;
//
//adjList[7].resize(2);
//adjList[7][0] = 5;
//adjList[7][1] = 8;
//
//adjList[8].resize(1);
//adjList[8][0] = 5;
//
//graph::StronglyConnectedComponentR<ID> sccr;
//sccr.solve(adjList);
//
//graph::StronglyConnectedComponentNR<ID> sccnr;
//sccnr.solve(adjList);
//
//if (sccr.sccs.size() == sccnr.sccs.size()) {
//    for (size_t i = 0; i < sccr.sccs.size(); ++i) {
//        if (sccr.sccs[i].size() == sccnr.sccs[i].size()) {
//            for (size_t j = 0; j < sccr.sccs[i].size(); ++j) {
//                if (sccr.sccs[i][j] != sccnr.sccs[i][j]) {
//                    cerr << "not match!" << endl;
//                }
//            }
//        } else {
//            cerr << "not match!" << endl;
//        }
//    }
//} else {
//    cerr << "not match!" << endl;
//}
