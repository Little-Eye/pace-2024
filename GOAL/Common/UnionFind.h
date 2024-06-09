////////////////////////////////
/// usage : 1.	union-find (a.k.a. disjoint set) data structure.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_UNION_FIND_H
#define CN_HUST_GOAL_COMMON_UNION_FIND_H


#include "GOAL/Typedef.h"


namespace goal {

struct UnionFind {
	Vec<ID> parents; // `parents[i]` is the parent of the item `i`.

	void init(ID itemNum) {
		parents.resize(itemNum);
		for (ID i = 0; i < itemNum; ++i) { parents[i] = i; }
	}

	ID parent(ID item) { // OPT[szx][8]: avoid stack overflow.
		return (item != parents[item]) ? (parents[item] = parent(parents[item])) : item;
	}

	bool join(ID item0, ID item1) { // return true if two sets are merged.
		ID parent0 = parent(item0);
		ID parent1 = parent(item1);
		if (parent0 == parent1) { return false; }
		parents[parent1] = parent0;
		return true;
	}
};

}


#endif // CN_HUST_GOAL_COMMON_UNION_FIND_H
