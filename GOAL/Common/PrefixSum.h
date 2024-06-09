////////////////////////////////
/// usage : 1.	querying the prefix sum and modifying item in O(log n) time.
/// 
/// note  : 1.	https://zhuanlan.zhihu.com/p/25185969
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_PREFIX_SUM_H
#define CN_HUST_GOAL_COMMON_PREFIX_SUM_H


#include "GOAL/Typedef.h"
#include "GOAL/Common/Math.h"


namespace goal {

// the index begins from 1.
class PrefixSum {
public:
    PrefixSum(ID count) : len(count + 1), mid(math::mostSignificantBit(count - 1)), segmentSums(len, 0) {}

    void increase(ID index, ID num) {
        do {
            segmentSums[index] += num;
        } while ((index += math::leastSignificantBit(index)) < len);
    }
    void increase_safe(ID index, ID num) { // make sure the index is in valid range.
        if (index < 1) { throw "invalid index"; }
        for (; index < len; index += math::leastSignificantBit(index)) {
            segmentSums[index] += num;
        }
    }

    ID accumulate(ID index) const {
        ID sum = segmentSums[index];
        while ((index -= math::leastSignificantBit(index)) > 0) {
            sum += segmentSums[index];
        }
        return sum;
    }

    // require each number at each index is non-negative.
    ID locate(ID sum) const {
        ID r = 0;
        for (ID mask = mid; mask > 0; mask >>= 1) {
            if (segmentSums[r | mask] < sum) {
                sum -= segmentSums[r |= mask];
            }
        }
        return r;
    }

protected:
    ID len;
    ID mid;
    Vec<ID> segmentSums;
};

// the index begins from 1 and the length is always the power of 2.
class PrefixSumPow2 : public PrefixSum {
public:
    PrefixSumPow2(ID exponent) : PrefixSum((1 << exponent) + 1) {}

    void increase(ID index, ID num) {
        segmentSums[index] += num;
        for (ID lb; (lb = math::leastSignificantBit(index)) < index;) {
            segmentSums[index += lb] += num;
        }
        for (; index < len; index <<= 1) {
            segmentSums[index] += num;
        }
    }
};

}


#endif // CN_HUST_GOAL_COMMON_PREFIX_SUM_H
