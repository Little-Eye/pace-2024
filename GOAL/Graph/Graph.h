////////////////////////////////
/// usage : 1.	topological or geometrical graph representations.
/// 
/// note  : 1.	it seems that `sizeof(WeightedAdjNode) == sizeof(ID) * 2`, which is inefficient and cause bugs in serialization.
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_GRAPH_H
#define CN_HUST_GOAL_COMMON_GRAPH_H


#include <initializer_list>

#include "GOAL/Typedef.h"
#include "GOAL/Common/Arr.h"


namespace goal {
namespace graph {

using DefaultWeightType = int;
using DefaultCapacityType = int;
using DefaultCoordType = int;


template<typename Weight>
struct AttrWeight {
	Weight weight;
};

template<typename Capacity>
struct AttrCapacity {
	Capacity capacity;
};

template<typename Weight, typename Capacity = DefaultCapacityType>
struct AttrWeightCapacity {
	Weight weight;
	Capacity capacity;
};

template<typename NodeId = ID, typename Attr = Void>
struct AdjNode : public Attr { // adjacent node.
	NodeId dst;
};

template<typename NodeId = ID, typename Weight = DefaultWeightType>
using WeightedAdjNode = AdjNode<NodeId, AttrWeight<Weight>>;

template<typename NodeId = ID, typename Capacity = DefaultCapacityType>
using CapacitatedAdjNode = AdjNode<NodeId, AttrCapacity<Capacity>>;

template<typename NodeId = ID, typename Weight = DefaultWeightType, typename Capacity = DefaultCapacityType>
using WeightedCapacitatedAdjNode = AdjNode<NodeId, AttrWeightCapacity<Weight, Capacity>>;


template<typename NodeId = ID, typename Attr = Void>
struct Arc : public Attr {
	NodeId src;
	NodeId dst;
};

template<typename NodeId = ID, typename Weight = DefaultWeightType>
using WeightedArc = Arc<NodeId, AttrWeight<Weight>>;

template<typename NodeId = ID, typename Capacity = DefaultCapacityType>
using CapacitatedArc = Arc<NodeId, AttrCapacity<Capacity>>;

template<typename NodeId = ID, typename Weight = DefaultWeightType, typename Capacity = DefaultCapacityType>
using WeightedCapacitatedArc = Arc<NodeId, AttrWeightCapacity<Weight, Capacity>>;


#pragma warning(push)
#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
// geometric graph elements.
template<typename Coord = DefaultCoordType>
using Coord2D = Arr<Coord, 2>;

template<typename Coord = DefaultCoordType>
using Coord3D = Arr<Coord, 3>;
#pragma warning(pop)

template<typename Coord = DefaultCoordType>
using Coords = Vec<Coord>;

template<typename Coord = DefaultCoordType>
inline Coord manhattanEdgeDistance(const Coord2D<Coord>& a, const Coord2D<Coord>& b) {
	return std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]);
}
template<typename Coord = DefaultCoordType>
inline Coord manhattanEdgeDistance(const Coord3D<Coord>& a, const Coord3D<Coord>& b) {
	return std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) + std::abs(a[2] - b[2]);
}
template<typename Coord = DefaultCoordType>
inline Coord manhattanEdgeDistance(const Coords<Coord>& a, const Coords<Coord>& b) {
	Coord d = 0;
	for (size_t c = 0; c < a.size(); ++c) { d += std::abs(a[c] - b[c]); }
	return d;
}
template<typename Coord = DefaultCoordType>
inline Coord manhattanNodeDistance(const Coord2D<Coord>& a, const Coord2D<Coord>& b) {
	return manhattanEdgeDistance(a, b) + 1;
}
template<typename Coord = DefaultCoordType>
inline Coord manhattanNodeDistance(const Coord3D<Coord>& a, const Coord3D<Coord>& b) {
	return manhattanEdgeDistance(a, b) + 1;
}
template<typename Coord = DefaultCoordType>
inline Coord manhattanNodeDistance(const Coords<Coord>& a, const Coords<Coord>& b) {
	return manhattanEdgeDistance(a, b) + 1;
}

template<typename Coord = DefaultCoordType>
inline Coord edgeVolume(const Coord2D<Coord>& a, const Coord2D<Coord>& b) {
	return std::abs(a[0] - b[0]) * std::abs(a[1] - b[1]);
}
template<typename Coord = DefaultCoordType>
inline Coord edgeVolume(const Coord3D<Coord>& a, const Coord3D<Coord>& b) {
	return std::abs(a[0] - b[0]) * std::abs(a[1] - b[1]) * std::abs(a[2] - b[2]);
}
template<typename Coord = DefaultCoordType>
inline Coord edgeVolume(const Coords<Coord>& a, const Coords<Coord>& b) {
	Coord d = 0;
	for (size_t c = 0; c < a.size(); ++c) { d *= std::abs(a[c] - b[c]); }
	return d;
}
template<typename Coord = DefaultCoordType>
inline Coord nodeVolume(const Coord2D<Coord>& a, const Coord2D<Coord>& b) {
	return (std::abs(a[0] - b[0]) + 1) * (std::abs(a[1] - b[1]) + 1);
}
template<typename Coord = DefaultCoordType>
inline Coord nodeVolume(const Coord3D<Coord>& a, const Coord3D<Coord>& b) {
	return (std::abs(a[0] - b[0]) + 1) * (std::abs(a[1] - b[1]) + 1) * (std::abs(a[2] - b[2]) + 1);
}
template<typename Coord = DefaultCoordType>
inline Coord nodeVolume(const Coords<Coord>& a, const Coords<Coord>& b) {
	Coord d = 0;
	for (size_t c = 0; c < a.size(); ++c) { d *= (std::abs(a[c] - b[c]) + 1); }
	return d;
}


// composite graph elements.
template<typename Weight = DefaultWeightType>
struct Tour {
	Weight distance;
	Vec<ID> nodes; // `nodes[i]` is the `i`_th node visited in the tour.
};
template<typename Weight = DefaultWeightType>
using Path = Tour<Weight>;


// graph representations.
template<typename T = DefaultWeightType, typename Index = int> // `T` can be weight, capacity, both of them or any other adjacency info.
using AdjMat = Mat<T, Index>; // adjacency matrix. `adjMat[n][m]` is the `T` on the arc from node `n` to node `m`.

template<typename AdjNodeType = WeightedAdjNode<>>
using AdjNodes = Vec<AdjNodeType>;

template<typename AdjNodeType = WeightedAdjNode<>>
using AdjList = Vec<AdjNodes<AdjNodeType>>; // adjacency list. `adjList[n][i]` is the `i`_th adjacent node to node `n`.

template<typename ArcType = WeightedArc<>>
using ArcList = Vec<ArcType>;

template<typename CoordType = Coord2D<>>
using CoordList = Vec<CoordType>;

}
}


#endif // CN_HUST_GOAL_COMMON_GRAPH_H
