////////////////////////////////
/// usage : 1.	depth first search framework.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_GRAPH_DEPTH_FIRST_SEARCH_H
#define CN_HUST_GOAL_GRAPH_DEPTH_FIRST_SEARCH_H


#include <functional>
#include <type_traits>

#include "GOAL/Typedef.h"
#include "GOAL/Common/LoopQueue.h"
#include "GOAL/Graph/Graph.h"


namespace goal {
namespace graph {
namespace dfs {

template<typename NodeId = ID>
using VisitNodeOnPush = Func<bool(NodeId to, NodeId from)>; // return true if the search should stop at this node (`dst`).

template<typename NodeId = ID>
using VisitNodeOnPop = Func<bool(NodeId to, NodeId from)>; // return true if the search should stop at this node (`dst`).

template<typename NodeId = ID>
using NextNeighbor = Func<NodeId(NodeId from)>; // return `InvalidId` if no more neighbors.

template<typename NodeId = ID>
using IsVisitedNode = Func<bool(NodeId to, NodeId from)>; // return true if the node has already been visited.

template<typename NodeId = ID>
void depthFirstSearch(Vec<NodeId>& s,
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	NextNeighbor<NodeId> nextNeighbor,
	IsVisitedNode<NodeId> isVisitedNode) {
	for (;;) {
		ID from = s.back();
		ID to = nextNeighbor(from);
		if (to != InvalidId) {
			if (isVisitedNode(to, from)) { continue; }
			if (visitNodeOnPush(to, from)) { return; }
			s.push_back(to);
		} else {
			s.pop_back();
			if (s.empty() || visitNodeOnPop(from, s.back())) { return; }
		}
	}
}

template<typename NodeId = ID>
void depthFirstSearch(Vec<NodeId>& s, NodeId src, // the source node of the search will be skipped (`visitNodeOnPush(?, src)` will not be called).
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	NextNeighbor<NodeId> nextNeighbor,
	IsVisitedNode<NodeId> isVisitedNode) {
	s.push_back(src);
	depthFirstSearch(s, visitNodeOnPush, visitNodeOnPop, nextNeighbor, isVisitedNode);
}

template<typename NodeId = ID>
void depthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search will be skipped (`visitNodeOnPush(?, src)` will not be called).
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	NextNeighbor<NodeId> nextNeighbor,
	IsVisitedNode<NodeId> isVisitedNode) {
	Vec<NodeId> s;
	s.reserve(nodeNum);
	depthFirstSearch<NodeId>(s, src, visitNodeOnPush, visitNodeOnPop, nextNeighbor, isVisitedNode);
}

template<typename NodeId = ID>
void depthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search will be skipped (`visitNodeOnPush(?, src)` will not be called).
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	NextNeighbor<NodeId> nextNeighbor) {
	if (std::is_integral<NodeId>::value) {
		Vec<bool> isVisited(nodeNum, false);
		isVisited[src] = true;
		depthFirstSearch<NodeId>(nodeNum, src,
			[&](NodeId to, NodeId from) { isVisited[to] = true; return visitNodeOnPush(to, from); },
			visitNodeOnPop, nextNeighbor,
			[&](NodeId to, NodeId) { return isVisited[to]; });
	} else {
		Set<NodeId> visitedNodes;
		visitedNodes.insert(src);
		depthFirstSearch<NodeId>(nodeNum, src,
			[&](NodeId to, NodeId from) { visitedNodes.insert(to); return visitNodeOnPush(to, from); },
			visitNodeOnPop, nextNeighbor,
			[&](NodeId to, NodeId) { return visitedNodes.find(to) != visitedNodes.end(); });
	}
}

template<typename AdjacentNode, typename NodeId = ID>
void depthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search (`visitNodeOnPush(?, src)` will not be called).
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	AdjList<AdjacentNode>& adjList,
	IsVisitedNode<NodeId> isVisitedNode) {
	Vec<ID> neighborIters(nodeNum, 0);
	depthFirstSearch<NodeId>(nodeNum, src, visitNodeOnPush, visitNodeOnPop,
		[&](NodeId from) {
			if (neighborIters[from] >= adjList[from].size()) { return InvalidId; }
			return adjList[from][neighborIters[from]++].dst;
		}, isVisitedNode);
}

template<typename AdjacentNode, typename NodeId = ID>
void depthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search (`visitNodeOnPush(?, src)` will not be called).
	VisitNodeOnPush<NodeId> visitNodeOnPush,
	VisitNodeOnPop<NodeId> visitNodeOnPop,
	AdjList<AdjacentNode>& adjList) {
	Vec<ID> neighborIters(nodeNum, 0);
	depthFirstSearch<NodeId>(nodeNum, src, visitNodeOnPush, visitNodeOnPop,
		[&](NodeId from) {
			if (neighborIters[from] >= adjList[from].size()) { return InvalidId; }
			return adjList[from][neighborIters[from]++].dst;
		});
}

}
}
}


#endif // CN_HUST_GOAL_GRAPH_DEPTH_FIRST_SEARCH_H

// sample code.
//using AdjNode = graph::AdjNode<ID>;
//ID n = 9;
//graph::AdjList<AdjNode> adjList(n);
//adjList[0].resize(1);
//adjList[0][0].dst = 1;
//
//adjList[1].resize(2);
//adjList[1][0].dst = 2;
//adjList[1][1].dst = 5;
//
//adjList[2].resize(2);
//adjList[2][0].dst = 3;
//adjList[2][1].dst = 4;
//
//adjList[3].resize(1);
//adjList[3][0].dst = 1;
//
//adjList[5].resize(2);
//adjList[5][0].dst = 4;
//adjList[5][1].dst = 6;
//
//adjList[6].resize(1);
//adjList[6][0].dst = 7;
//
//adjList[7].resize(2);
//adjList[7][0].dst = 5;
//adjList[7][1].dst = 8;
//
//adjList[8].resize(1);
//adjList[8][0].dst = 5;
//
//graph::dfs::depthFirstSearch<AdjNode, ID>(n, 0,
//	[](ID t, ID f) { cout << "i" << t << endl; return false; },
//	[](ID t) { cout << "o" << t << endl; return false; }, adjList);
