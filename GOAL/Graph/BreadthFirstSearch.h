////////////////////////////////
/// usage : 1.	breadth first search framework.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_GRAPH_BREADTH_FIRST_SEARCH_H
#define CN_HUST_GOAL_GRAPH_BREADTH_FIRST_SEARCH_H


#include <functional>
#include <type_traits>

#include "GOAL/Typedef.h"
#include "GOAL/Common/LoopQueue.h"
#include "GOAL/Graph/Graph.h"


namespace goal {
namespace graph {
namespace bfs {

template<typename NodeId = ID>
using VisitNode = Func<bool(NodeId to, NodeId from)>; // return true if the search should stop at this node (`dst`).

template<typename NodeId = ID>
using TraverseNeighbors = Func<bool(NodeId from, VisitNode<NodeId> onNeighbor)>; // return true immediately if `onNeighbor` returns true.

template<typename NodeId = ID>
using IsVisitedNode = Func<bool(NodeId)>; // return true if the node has already been visited.

template<typename NodeId = ID>
void breadthFirstSearch(LoopQueue<NodeId, NodeId>& q,
	VisitNode<NodeId> visitNode,
	TraverseNeighbors<NodeId> traverseNeighbors,
	IsVisitedNode<NodeId> isVisitedNode) {
	auto onNeighbor = [&](NodeId to, NodeId from) {
		if (isVisitedNode(to)) { return false; }
		if (visitNode(to, from)) { return true; }
		q.pushBack(to);
		return false;
	};
	for (; !q.empty(); q.popFront()) {
		if (traverseNeighbors(q.front(), onNeighbor)) { return; }
	}
}

template<typename NodeId = ID>
void breadthFirstSearch(LoopQueue<NodeId, NodeId>& q, NodeId src, // the source node of the search will be skipped (`visitNode(src)` will not be called).
	VisitNode<NodeId> visitNode,
	TraverseNeighbors<NodeId> traverseNeighbors,
	IsVisitedNode<NodeId> isVisitedNode) {
	q.pushBack(src);
	breadthFirstSearch(q, visitNode, traverseNeighbors, isVisitedNode);
}

template<typename NodeId = ID>
void breadthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search will be skipped (`visitNode(src)` will not be called).
	VisitNode<NodeId> visitNode,
	TraverseNeighbors<NodeId> traverseNeighbors,
	IsVisitedNode<NodeId> isVisitedNode) {
	LoopQueue<NodeId, NodeId> lq(nodeNum);
	breadthFirstSearch<NodeId>(lq, src, visitNode, traverseNeighbors, isVisitedNode);
}

template<typename NodeId = ID>
void breadthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search will be skipped (`visitNode(src)` will not be called).
	VisitNode<NodeId> visitNode,
	TraverseNeighbors< NodeId> traverseNeighbors) {
	if (std::is_integral<NodeId>::value) {
		Vec<bool> isVisited(nodeNum, false);
		isVisited[src] = true;
		breadthFirstSearch<NodeId, NodeId>(nodeNum, src,
			[&](NodeId to, NodeId from) { isVisited[to] = true; return visitNode(to, from); },
			traverseNeighbors,
			[&](NodeId to) { return isVisited[to]; });
	} else {
		Set<NodeId> visitedNodes;
		visitedNodes.insert(src);
		breadthFirstSearch<NodeId, NodeId>(nodeNum, src,
			[&](NodeId to, NodeId from) { visitedNodes.insert(to); return visitNode(to, from); },
			traverseNeighbors,
			[&](NodeId to) { return visitedNodes.find(to) != visitedNodes.end(); });
	}
}

template<typename AdjacentNode, typename NodeId = ID>
void breadthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search (`visitNode(src)` will not be called).
	VisitNode<NodeId> visitNode,
	const AdjList<AdjacentNode>& adjList,
	IsVisitedNode<NodeId> isVisitedNode) {
	breadthFirstSearch<NodeId>(nodeNum, src, visitNode,
		[&](NodeId from, VisitNode<NodeId> onNeighbor) {
			for (auto to = adjList[from].begin(); to != adjList[from].end(); ++to) {
				if (onNeighbor(to->dst, from)) { return true; }
			} // FIX[szx][9]: make sure `AdjacentNode` inherits `graph::AdjNode`.
			return false;
		}, isVisitedNode);
}

template<typename AdjacentNode, typename NodeId = ID>
void breadthFirstSearch(NodeId nodeNum, NodeId src, // the source node of the search (`visitNode(src)` will not be called).
	VisitNode<NodeId> visitNode,
	const AdjList<AdjacentNode>& adjList) {
	breadthFirstSearch<NodeId>(nodeNum, src, visitNode,
		[&](NodeId from, VisitNode<NodeId> onNeighbor) {
			for (auto to = adjList[from].begin(); to != adjList[from].end(); ++to) {
				if (onNeighbor(to->dst, from)) { return true; }
			} // FIX[szx][9]: make sure `AdjacentNode` inherits `graph::AdjNode`.
			return false;
		});
}

}
}
}


#endif // CN_HUST_GOAL_GRAPH_BREADTH_FIRST_SEARCH_H
