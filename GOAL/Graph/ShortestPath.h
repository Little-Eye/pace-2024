////////////////////////////////
/// usage : 1.	shortest (simple) path algorithms.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_GRAPH_SHORTEST_PATH_H
#define CN_HUST_GOAL_GRAPH_SHORTEST_PATH_H


#include "GOAL/Typedef.h"
#include "GOAL/Common/LoopQueue.h"
#include "GOAL/Graph/Graph.h"


namespace goal {
namespace graph {

namespace floyd {

// find shortest paths between each pair of nodes.
// src -> nextMedian[src, dst] ->->->->-> dst
template<typename Weight, typename NodeId = ID>
void findAllPairsPaths(Mat<Weight, NodeId>& adjMat, Mat<NodeId, NodeId>& nextMedian) {
	NodeId nodeNum = adjMat.size1();

	nextMedian.init(nodeNum, nodeNum);
	auto iter = nextMedian.begin();
	for (NodeId src = 0; src < nodeNum; ++src) {
		for (NodeId dst = 0; dst < nodeNum; ++dst) { *(iter++) = dst; }
	}

	for (NodeId mid = 0; mid < nodeNum; ++mid) {
		for (NodeId src = 0; src < nodeNum; ++src) {
			for (NodeId dst = 0; dst < nodeNum; ++dst) {
				Weight w = adjMat.at(src, mid) + adjMat.at(mid, dst);
				if (w < adjMat.at(src, dst)) {
					adjMat.at(src, dst) = w;
					nextMedian.at(src, dst) = nextMedian.at(src, mid);
				}
			}
		}
	}
}

template<typename Weight, typename NodeId = ID>
void findAllPairsPaths(Mat<Weight, NodeId>& adjMat) {
	NodeId nodeNum = adjMat.size1();

	for (NodeId mid = 0; mid < nodeNum; ++mid) {
		for (NodeId src = 0; src < nodeNum; ++src) {
			for (NodeId dst = 0; dst < nodeNum; ++dst) {
				Weight w = adjMat.at(src, mid) + adjMat.at(mid, dst);
				if (w < adjMat.at(src, dst)) { adjMat.at(src, dst) = w; }
			}
		}
	}
}

template<typename NodeId = ID>
Vec<NodeId> getPath(const Mat<NodeId, NodeId>& nextMedian, ID src, ID dst) {
	Vec<NodeId> path;
	path.reserve(nextMedian.size1());
	while (src != dst) { path.push_back(src = nextMedian.at(src, dst)); }
	return path;
}

}

namespace bellman {

template<typename Weight, typename NodeId = ID>
using GetCost = Func<Weight(NodeId to)>;

template<typename Weight, typename NodeId = ID>
using SetCost = Func<void(NodeId to, Weight cost)>;

template<typename NodeId = ID>
using SetPrevNode = Func<void(NodeId nextNode, NodeId thisNode)>;

template<typename Weight, typename NodeId = ID>
void findSingleSourcePaths(const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, NodeId src,
	GetCost<Weight, NodeId> getCost, SetCost<Weight, NodeId> setCost, SetPrevNode<NodeId> setPrevNode) {
	NodeId nodeNum = adjList.size();

	LoopQueue<NodeId, NodeId> nodesToRelax(nodeNum);
	Vec<bool> inQueue(nodeNum, false);

	setCost(src, 0);
	nodesToRelax.pushBack(src);
	inQueue[src] = true;

	while (!nodesToRelax.empty()) {
		NodeId node = nodesToRelax.front(); // OPTIMIZE[szx][8]: Large Label Last optimization.

		inQueue[node] = false;
		nodesToRelax.popFront();

		const AdjNodes<WeightedAdjNode<NodeId, Weight>>& adjNodes(adjList[node]);
		for (auto adjNode = adjNodes.begin(); adjNode != adjNodes.end(); ++adjNode) {
			Weight newDist = getCost(node) + adjNode->weight;
			if (newDist < getCost(adjNode->dst)) {
				setCost(adjNode->dst, newDist);
				setPrevNode(adjNode->dst, node);

				if (!inQueue[adjNode->dst]) {
					inQueue[adjNode->dst] = true;
					(getCost(adjNode->dst) < getCost(nodesToRelax.front())) // Small Label First optimization.
						? nodesToRelax.pushFront(adjNode->dst)
						: nodesToRelax.pushBack(adjNode->dst);
				}
			}
		}
	}
}

template<typename Weight, typename NodeId = ID, Weight MaxArcWeight = SafeLimit<Weight>::max>
void findSingleSourcePaths(const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, NodeId src, Vec<Weight>& minCost, Vec<NodeId>& prevNode) {
	NodeId nodeNum = adjList.size();

	minCost.clear();
	minCost.resize(nodeNum, MaxArcWeight);

	prevNode.clear();
	prevNode.resize(nodeNum, InvalidId);

	findSingleSourcePaths<Weight, NodeId>(adjList, src,
		[&minCost](NodeId node) { return minCost[node]; },
		[&minCost](NodeId node, Weight cost) { minCost[node] = cost; },
		[&prevNode](NodeId nextNode, NodeId thisNode) { prevNode[nextNode] = thisNode; });
}

template<typename Weight, typename NodeId = ID, Weight MaxArcWeight = SafeLimit<Weight>::max>
void findSingleSourcePaths(const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, NodeId src, Vec<Weight>& minCost) {
	NodeId nodeNum = adjList.size();

	minCost.clear();
	minCost.resize(nodeNum, MaxArcWeight);

	findSingleSourcePaths<Weight, NodeId>(adjList, src,
		[&minCost](NodeId node) { return minCost[node]; },
		[&minCost](NodeId node, Weight cost) { minCost[node] = cost; },
		[](NodeId, NodeId) {});
}

template<typename Weight, typename NodeId = ID>
void findAllPairsPaths(const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, Mat<Weight, NodeId>& adjMat, Mat<NodeId, NodeId>& prevNode) {
	NodeId nodeNum = adjList.size();

	adjMat.init(nodeNum, nodeNum);
	adjMat.reset(ArrResetOption::SafeMaxInt);

	prevNode.clear();
	prevNode.reset(nodeNum, InvalidId);

	for (NodeId src = 0; src < nodeNum; ++src) {
		typename Mat<Weight, NodeId>::Iterator minCost = adjMat[src];
		findSingleSourcePaths<Weight, NodeId>(adjList, src,
			[&minCost](NodeId node) { return minCost[node]; },
			[&minCost](NodeId node, Weight cost) { minCost[node] = cost; },
			[&prevNode](NodeId nextNode, NodeId thisNode) { prevNode[nextNode] = thisNode; });
	}
}

template<typename Weight, typename NodeId = ID>
void findAllPairsPaths(const AdjList<WeightedAdjNode<NodeId, Weight>>& adjList, Mat<Weight, NodeId>& adjMat) {
	NodeId nodeNum = adjList.size();

	adjMat.init(nodeNum, nodeNum);
	adjMat.reset(ArrResetOption::SafeMaxInt);

	for (NodeId src = 0; src < nodeNum; ++src) {
		typename Mat<Weight, NodeId>::Iterator minCost = adjMat[src];
		findSingleSourcePaths<Weight, NodeId>(adjList, src,
			[&minCost](NodeId node) { return minCost[node]; },
			[&minCost](NodeId node, Weight cost) { minCost[node] = cost; },
			[](NodeId, NodeId) {});
	}
}

template<typename NodeId = ID>
Vec<NodeId> getReversePath(const Vec<NodeId>& prevNode, ID dst, bool noSrc = true, bool noDst = false) {
	Vec<NodeId> path;
	path.reserve(prevNode.size());
	if (noDst) { dst = prevNode[dst]; }
	for (; dst != InvalidId; dst = prevNode[dst]) { path.push_back(dst); }
	if (noSrc) { path.pop_back(); }
	return path;
}

template<typename NodeId = ID>
Vec<NodeId> getPath(const Vec<NodeId>& prevNode, ID dst, bool noSrc = true, bool noDst = false) {
	Vec<NodeId> path;
	getReversePath(path, dst, noSrc, noDst);
	std::reverse(path.begin(), path.end());
	return path;
}

}

// EXT[szx][9]: add more shortest path algorihtm (https://github.com/Zhouxing-Su/CppUtilibs/blob/master/algorithm/Graph.h#L303).

}
}


#endif // CN_HUST_GOAL_GRAPH_SHORTEST_PATH_H
