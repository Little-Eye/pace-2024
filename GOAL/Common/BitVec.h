////////////////////////////////
/// usage:	1.	variable length bit vector.
/// 
/// note:	1.	
////////////////////////////////

#ifndef CN_HUST_HAFV_UTIL_BIT_VEC_H
#define CN_HUST_HAFV_UTIL_BIT_VEC_H


#include <algorithm>
#include <cstdint>
#include <climits>

#include "GOAL/Typedef.h"


namespace goal {

class BitVec {
public:
    using Index = size_t;


    static constexpr Index BitPerByte = CHAR_BIT;
    static constexpr Index BytePerBlock = sizeof(int);
    static constexpr Index BitPerBlock = BytePerBlock * BitPerByte;


    BitVec() {}
    BitVec(Index size) : bits(calcBlockNum(size)), count(size) {}
    BitVec(Index size, bool initValue) : BitVec(size) {
        std::fill(bits.begin(), bits.end(), (initValue ? -1 : 0));
    }

    void set(Index i) { bits[blockIndex(i)] |= (1 << blockOffset(i)); }
    void reset(Index i) { bits[blockIndex(i)] &= (~(1 << blockOffset(i))); }

    bool operator[](Index i) const { return (bits[blockIndex(i)] & (~(1 << blockOffset(i)))) != 0; }

    char* data() { return toBytes(bits.data()); }

    static Index calcBlockNum(Index size) { return (size + BitPerBlock - 1) / BitPerBlock; }
    static Index calcByteNum(Index size) { return (size + BitPerBlock - 1) / BitPerByte; }

protected:
    // which block.
    Index blockIndex(Index i) const { return i / BitPerBlock; }
    // which bit in block.
    Index blockOffset(Index i) const { return i % BitPerBlock; }


    Vec<int> bits;
    Index count;
};

}


#endif // CN_HUST_HAFV_UTIL_BIT_VEC_H
