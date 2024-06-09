////////////////////////////////
/// usage : 1.	hash-based element existence checker.
/// 
/// note  : 1.	false positive results may be obtained.
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMOM_BLOOM_FILTER
#define CN_HUST_GOAL_COMMOM_BLOOM_FILTER


#include "GOAL/Typedef.h"
#include "GOAL/Common/Arr.h"


namespace goal {

struct Hash {
    static constexpr int MersennePrime = 2147483647; // the max Mersenne prime within the range of int.

    // accumulator.
    template<typename HashVal = unsigned>
    static HashVal addAcc(HashVal &r, HashVal n) { return r += n; }

    template<typename HashVal = unsigned>
    static HashVal xorAcc(HashVal &r, HashVal n) { return r ^= n; }

    template<typename HashVal = unsigned>
    static HashVal multAcc(HashVal &r, HashVal n) { return r *= n; }

    // transformer.
    template<typename HashVal = unsigned>
    static HashVal asIs(HashVal n) { return n; }

    template<typename HashVal = unsigned>
    static HashVal twist(HashVal n) { return MersennePrime * n; }

    template<typename HashVal = unsigned>
    static HashVal square(HashVal n) { return n * n; }

    template<typename HashVal = unsigned>
    static HashVal cube(HashVal n) { return n * n * n; }

    template<typename HashVal = unsigned>
    static HashVal quad(HashVal n) { return square(n) * square(n); }

    // calculator.
    template<typename HashVal = unsigned, typename Accumulator, typename Transformer, typename HashValList = Vec<HashVal>>
    static HashVal calc(const HashValList &items, Accumulator accumulate, Transformer transform) {
        HashVal r = transform(items[0]);
        for (auto i = items.begin() + 1; i != items.end(); ++i) { accumulate(r, transform(*i)); }
        return r;
    }
    template<typename HashVal = unsigned, typename Accumulator, typename Transformer, typename HashValList = Vec<HashVal>>
    static HashVal calcCoef(const HashValList &items, Accumulator accumulate, Transformer transform) {
        HashVal c = 1;
        HashVal r = transform(items[0]);
        for (auto i = items.begin() + 1; i != items.end(); ++i) { accumulate(r, transform(*i * (++c))); }
        return r;
    }
    template<typename HashVal = unsigned, typename Accumulator, typename Transformer, typename HashValList = Vec<HashVal>>
    static HashVal calcDiff(const HashValList &items, Accumulator accumulate, Transformer transform) {
        HashVal prev = items[1];
        HashVal r = transform(prev - items[0]);
        for (auto i = items.begin() + 2; i != items.end(); prev = *i, ++i) { accumulate(r, transform(*i - prev)); }
        return r;
    }
    template<typename HashVal = unsigned, typename Accumulator, typename Transformer, typename HashValList = Vec<HashVal>>
    static HashVal calcCoefDiff(const HashValList &items, Accumulator accumulate, Transformer transform) {
        HashVal c = 1;
        HashVal prev = items[1];
        HashVal r = transform(prev - items[0]);
        for (auto i = items.begin() + 2; i != items.end(); prev = *i, ++i) { accumulate(r, transform((*i - prev) * (++c))); }
        return r;
    }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal directSum(const HashValList &items) { return calc<HashVal>(items, addAcc<HashVal>, asIs<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal twistedSum(const HashValList &items) { return calc<HashVal>(items, addAcc<HashVal>, twist<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal squareSum(const HashValList &items) { return calc<HashVal>(items, addAcc<HashVal>, square<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal cubeSum(const HashValList &items) { return calc<HashVal>(items, addAcc<HashVal>, cube<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal quadSum(const HashValList &items) { return calc<HashVal>(items, addAcc<HashVal>, quad<HashVal>); }

    // [permutation][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefSum(const HashValList &items) { return calcCoef<HashVal>(items, addAcc<HashVal>, asIs<HashVal>); }

    // [permutation][incremental]
    //template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    //static HashVal diffSum(const HashValList &items) { return calcDiff<HashVal>(items, addAcc<HashVal>, asIs<HashVal>); } // meaningless (becomes back() - front()).

    // [permutation][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefDiffSum(const HashValList &items) { return calcCoefDiff<HashVal>(items, addAcc<HashVal>, asIs<HashVal>); }


    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal directXor(const HashValList &items) { return calc<HashVal>(items, xorAcc<HashVal>, asIs<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal twistedXor(const HashValList &items) { return calc<HashVal>(items, xorAcc<HashVal>, twist<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal squareXor(const HashValList &items) { return calc<HashVal>(items, xorAcc<HashVal>, square<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal cubeXor(const HashValList &items) { return calc<HashVal>(items, xorAcc<HashVal>, cube<HashVal>); }

    // [combination][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal quadXor(const HashValList &items) { return calc<HashVal>(items, xorAcc<HashVal>, quad<HashVal>); }

    // [permutation][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefXor(const HashValList &items) { return calcCoef<HashVal>(items, xorAcc<HashVal>, asIs<HashVal>); }

    // [permutation][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal diffXor(const HashValList &items) { return calcDiff<HashVal>(items, xorAcc<HashVal>, asIs<HashVal>); }

    // [permutation][incremental]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefDiffXor(const HashValList &items) { return calcCoefDiff<HashVal>(items, xorAcc<HashVal>, asIs<HashVal>); }


    // [combination]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal directMult(const HashValList &items) { return calc<HashVal>(items, multAcc<HashVal>, asIs<HashVal>); }

    // [combination]
    //template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    //static HashVal twistedMult(const HashValList &items) { return calc<HashVal>(items, multAcc<HashVal>, twist<HashVal>); } // meaningless (multiply a constant to each item).

    // [combination]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal squareMult(const HashValList &items) { return calc<HashVal>(items, multAcc<HashVal>, square<HashVal>); }

    // [combination]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal cubeMult(const HashValList &items) { return calc<HashVal>(items, multAcc<HashVal>, cube<HashVal>); }

    // [combination]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal quadMult(const HashValList &items) { return calc<HashVal>(items, multAcc<HashVal>, quad<HashVal>); }

    // [permutation]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefMult(const HashValList &items) { return calcCoef<HashVal>(items, multAcc<HashVal>, asIs<HashVal>); }

    // [permutation]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal diffMult(const HashValList &items) { return calcDiff<HashVal>(items, multAcc<HashVal>, asIs<HashVal>); }

    // [permutation]
    template<typename HashVal = unsigned, typename HashValList = Vec<HashVal>>
    static HashVal coefDiffMult(const HashValList &items) { return calcCoefDiff<HashVal>(items, multAcc<HashVal>, asIs<HashVal>); }
};


template<typename HashValue = unsigned>
class BloomFilter {
public:
    using HashVal = HashValue;


    BloomFilter(HashVal hashFuncNum, HashVal modulo) : flags(hashFuncNum, Vec<bool>(modulo)) {
        reset();
    }


    void reset() {
        for (auto f = flags.begin(); f != flags.end(); ++f) {
            std::fill(f->begin(), f->end(), false);
        }
        #if SZX_DEBUG
        queryCount = 0;
        hitCount = 0;
        insertCount = 0;
        nonZeroCount = 0;
        #endif // SZX_DEBUG
    }

    HashVal hashFuncNum() const { return sCast<HashVal>(flags.size()); }
    HashVal modulo() const { return sCast<HashVal>(flags[0].size()); }

    template<typename HashValList = Vec<HashVal>>
    bool contain(const HashValList &hashValues) const {
        #if SZX_DEBUG
        ++queryCount;
        #endif // SZX_DEBUG
        for (HashVal h = 0; h < hashFuncNum(); ++h) {
            if (!flags[h][hashValues[h]]) { return false; }
        }
        #if SZX_DEBUG
        ++hitCount;
        #endif // SZX_DEBUG
        return true;
    }

    template<typename HashValList = Vec<HashVal>>
    void insert(const HashValList &hashValues) {
        #if SZX_DEBUG
        ++insertCount;
        #endif // SZX_DEBUG
        for (HashVal h = 0; h < hashFuncNum(); ++h) {
            #if SZX_DEBUG
            if (!flags[h][hashValues[h]]) { ++nonZeroCount; }
            #endif // SZX_DEBUG
            flags[h][hashValues[h]] = true;
        }
    }

    #if SZX_DEBUG
    mutable long long queryCount;
    mutable long long hitCount;
    long long insertCount;
    long long nonZeroCount;
    #endif // SZX_DEBUG

protected:
    Vec<Vec<bool>> flags; // `flags[f][v] == true` means the items with hash value `v` by hash function `f` exists.
};

}


#endif // CN_HUST_GOAL_COMMOM_BLOOM_FILTER
