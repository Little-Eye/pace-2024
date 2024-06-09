////////////////////////////////
/// usage : 1.	fundamental parallel algorithms.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_PARALLEL_ALGORITHM_H
#define CN_HUST_GOAL_PARALLEL_ALGORITHM_H


#include "GOAL/Common/JobFlow.h"
#include "GOAL/Common/Math.h"
#include "GOAL/Common/Arr.h"


namespace goal {

// return true if any elements in `l` are true.
template<typename BoolList, typename IndexType = int>
bool fanAnd(JobFlow &jf, const BoolList &l) {
    bool r = true;
    jf.run([&](int dim, int idx) {
        for (IndexType i = idx; i < l.size(); i += dim) {
            if (!l[i]) { r = false; break; }
        }
    });
    return r;
}

// return true if any elements in `l` are true.
template<typename BoolList, typename IndexType = int>
bool fanOr(JobFlow &jf, const BoolList &l) {
    bool r = false;
    jf.run([&](int dim, int idx) {
        for (IndexType i = idx; i < l.size(); i += dim) {
            if (l[i]) { r = true; break; }
        }
    });
    return r;
}

template<typename Value, typename ValueList, typename IndexType = int>
Value fanMax4Serial(ValueList &l, IndexType size) {
    switch (size) {
    case 1: return l[0];
    case 2: return (std::max)(l[0], l[1]);
    case 3: return (std::max)(l[0], (std::max)(l[1], l[2]));
    case 4:
    default: return (std::max)((std::max)(l[0], l[1]), (std::max)(l[2], l[3]));
    }
}
// it may be better choice when `jf.workerNum()` is far less than `l.size()`.
template<typename Value, typename ValueList, typename IndexType = int>
Value fanMaxInPlace(JobFlow &jf, ValueList &l) {
    IndexType size = sCast<IndexType>(l.size());
    if (size <= 4) { return fanMax4Serial(l, size); }

    IndexType blockNum = (std::min)(jf.workerNum(), size / 2);
    jf.run([&](int dim, int idx) {
        if (idx >= blockNum) { return; }
        for (IndexType i = idx; (i += dim) < size;) {
            math::updateMax(l[idx], l[i]);
        }
    });
    return fanMaxInPlace(jf, maxValues);
}
// it may be better choice when `jf.workerNum()` is far less than `l.size()`.
template<typename Value, typename ValueList, typename IndexType = int>
Value fanMax(JobFlow &jf, const ValueList &l) {
    IndexType size = sCast<IndexType>(l.size());
    if (size <= 4) { return fanMax4Serial(l, size); }

    IndexType blockNum = (std::min)(jf.workerNum(), size / 2);
    ValueList maxValues(blockNum);
    jf.run([&](int dim, int idx) {
        if (idx >= blockNum) { return; }
        maxValues[idx] = l[idx];
        for (IndexType i = idx; (i += dim) < size;) {
            math::updateMax(maxValues[idx], l[i]);
        }
    });
    return fanMaxInPlace(jf, maxValues);
} // OPT[szx][9]: O(log log n) time O(n) workers implementation. O(1) time O(n^2) workers implementation.

}


#endif // CN_HUST_GOAL_PARALLEL_ALGORITHM_H
