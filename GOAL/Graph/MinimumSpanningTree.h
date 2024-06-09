////////////////////////////////
/// usage : 1.	minimum spanning tree algorithms.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_GRAPH_MINIMUM_SPANNING_TREE_H
#define CN_HUST_GOAL_GRAPH_MINIMUM_SPANNING_TREE_H


#include "GOAL/Typedef.h"
#include "GOAL/Graph/Graph.h"
#include "GOAL/Common/PriorityQueue.h"
#include "GOAL/Common/Math.h"


namespace goal {
namespace graph {
namespace mst {

// Prim's algorithm.
template<typename NodeId, typename Weight>
Weight minimumSpanningTree(Vec<NodeId>& prevNodes, const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, Weight maxWeight, NodeId root) {
	prevNodes.resize(adjList.size());

	List<Weight> prevNodesWeight(adjList.size());
	if (maxWeight <= SafeLimit<NodeId>::max) {
		prevNodesWeight.reset(ArrResetOption::SafeMaxInt);
	} else {
		std::fill(prevNodesWeight.begin(), prevNodesWeight.end(), maxWeight);
	}
	prevNodes[root] = root;
	prevNodesWeight[root] = 0;

	Weight w = 0;

	PriorityQueueByBucketL1<NodeId> edgePQ(maxWeight);
	//edgePQ.reserve(adjList.size() / 2);
	for (auto n = adjList[root].begin(); n != adjList[root].end(); ++n) {
		edgePQ.push({ root, n->dst }, n->weight);
		prevNodes[n->dst] = root;
		prevNodesWeight[n->dst] = n->weight;
	}
	for (NodeId restNodeNum = sCast<NodeId>(adjList.size()) - 1; !edgePQ.empty(); edgePQ.pop()) {
		NodeId nextNode = edgePQ.top(); // next node to span.
		if (prevNodesWeight[nextNode] < maxWeight) { continue; } // skip already spanned nodes.

		w += edgePQ.topPriority();
		if (--restNodeNum <= 0) { break; }

		for (auto n = adjList[nextNode].begin(); n != adjList[nextNode].end(); ++n) {
			if (math::updateMin(prevNodesWeight[n->dst], n->weight)) {
				prevNodes[n->dst] = nextNode;
				edgePQ.push(n->dst, n->weight);
			}
		}

	}

	return w;
}

}
}
}


#endif // CN_HUST_GOAL_GRAPH_MINIMUM_SPANNING_TREE_H
