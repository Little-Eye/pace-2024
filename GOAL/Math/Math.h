////////////////////////////////
/// usage : 1.	basic arithmetic and logic functions which prefers speed than accuracy.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_MATH_H
#define CN_HUST_GOAL_COMMON_MATH_H


#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>

#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"


namespace goal {

template<typename T>
struct SafeLimit { // no overflow occurs after single binary additive operation.
	static constexpr T max = sCast<T>(0x3F3F3F3F3F3F3F3Fll); //(std::numeric_limits<T>::max)() / 2 - 1;
};

template<typename T>
struct NumErr { static constexpr T v = 0; };
template<>
struct NumErr<float> { static constexpr float v = 1E-4f; };
template<>
struct NumErr<double> { static constexpr double v = 1E-8; };

template<typename IntegerType = int>
struct BitMask {
	static constexpr IntegerType DropLastBit = IntegerType(-1) ^ IntegerType(1); // least significant bit.
};

namespace math {

static constexpr Real Pi = 3.14159265358979323846;
static constexpr Real HalfPi = 1.57079632679489661923;
static constexpr Real Pi2 = 6.28318530717958647692;

#pragma region SetOperation
// setSymDiff(a, b) = |union(a, b) - intersetion(a, b)|.
// `orderedList1` and `orderedList2` are non-decreasing lists of numbers.
template<typename Iter, typename Size = int>
Size setSymDiffSize(Iter orderedList1Begin, Iter orderedList1End, Iter orderedList2Begin, Iter orderedList2End) {
	Size diff = 0;
	while (orderedList1Begin != orderedList1End) {
		if (orderedList2Begin == orderedList2End) {
			return diff + sCast<Size>(orderedList1End - orderedList1Begin);
		}

		//if (*i1 == *i2) { ++i1; ++i2; }
		//else if (*i1 < *i2) { ++diff; ++i1; }
		//else { ++diff; ++i2; }
		if (*orderedList1Begin < *orderedList2Begin) {
			++diff;
			++orderedList1Begin;
		} else {
			if (*orderedList2Begin < *orderedList1Begin) {
				++diff;
			} else {
				++orderedList1Begin;
			}
			++orderedList2Begin;
		}
	}
	return diff + sCast<Size>(orderedList2End - orderedList2Begin);
}

// `orderedList1` and `orderedList2` are non-decreasing lists of numbers.
template<typename List>
void setSymDiff(const List& orderedList1, const List& orderedList2, List& diff1, List& diff2) {
	bool shorterFirst = orderedList1.back() < orderedList2.back();
	const List& shorter(shorterFirst ? orderedList1 : orderedList2);
	const List& longer(!shorterFirst ? orderedList1 : orderedList2);
	List& shorterDiff(shorterFirst ? diff1 : diff2);
	List& longerDiff(!shorterFirst ? diff1 : diff2);
	auto l = longer.begin();
	for (auto s = shorter.begin(); s != shorter.end();) {
		if (*s == *l) {
			++s;
			++l;
		} else if (*l < *s) {
			longerDiff.push_back(*(l++));
		} else { // if (*s < *l)
			shorterDiff.push_back(*(s++));
		}
	}
	longerDiff.insert(longerDiff.end(), l, longer.end());
}
template<typename List>
void setSymDiff(const List& orderedList1, const List& orderedList2, List& diff1, List& diff2, List& common) {
	bool shorterFirst = orderedList1.back() < orderedList2.back();
	const List& shorter(shorterFirst ? orderedList1 : orderedList2);
	const List& longer(!shorterFirst ? orderedList1 : orderedList2);
	List& shorterDiff(shorterFirst ? diff1 : diff2);
	List& longerDiff(!shorterFirst ? diff1 : diff2);
	auto l = longer.begin();
	for (auto s = shorter.begin(); s != shorter.end();) {
		if (*s == *l) {
			++s;
			common.push_back(*(l++));
		} else if (*l < *s) {
			longerDiff.push_back(*(l++));
		} else { // if (*s < *l)
			shorterDiff.push_back(*(s++));
		}
	}
	longerDiff.insert(longerDiff.end(), l, longer.end());
}

// `common` and `orderedList` are non-decreasing lists of numbers.
template<typename List>
void setCommon(List& common, const List& orderedList) {
	while (!common.empty() && (common.back() > orderedList.back())) { common.pop_back(); }
	auto l = orderedList.begin();
	auto s = common.begin();
	for (; (s != common.end()) && (*s == *l); ++s, ++l) {}
	auto c = s;
	while (s != common.end()) {
		if (*s == *l) {
			*(c++) = *(s++);
			++l;
		} else if (*l < *s) {
			++l;
		} else { // if (*s < *l)
			++s;
		}
	}
	common.erase(c, common.end());
}

// `minuend` and `subtrahend` are non-decreasing lists of numbers.
// each item in `subtrahend` must also be in `minuend`.
template<typename List>
void setMinus(List& minuend, const List& subtrahend) {
	auto l = minuend.begin();
	auto s = subtrahend.begin();
	for (; (s != subtrahend.end()) && (*l < *s); ++l) {}
	auto m = l;
	while (s != subtrahend.end()) {
		if (*s == *l) {
			++s; ++l; // skip the common item.
		} else if (*l < *s) {
			*(m++) = *(l++); // save the dedicated item in `minuend`.
		} else { // if (*s < *l)
			++s; // skip the dedicated item in `subtrahend`.
		}
	}
	while (l < minuend.end()) { *(m++) = *(l++); }
	minuend.erase(m, minuend.end());
}
#pragma endregion SetOperation

#pragma region Compare
template<typename RealNum>
bool weakEqual(RealNum l, RealNum r) {
	RealNum d = l - r;
	return (-NumErr<RealNum>::v <= d) && (d <= NumErr<RealNum>::v);
}
template<typename RealNum>
bool weakLess(RealNum l, RealNum r) { // operator<=().
	return ((l - r) <= NumErr<RealNum>::v);
}
template<typename RealNum>
bool weakGreater(RealNum l, RealNum r) { // operator>=().
	return ((l - r) >= -NumErr<RealNum>::v);
}
template<typename RealNum>
bool strictLess(RealNum l, RealNum r) { // operator<().
	return ((l - r) < -NumErr<RealNum>::v);
}
template<typename RealNum>
bool strictGreater(RealNum l, RealNum r) { // operator>().
	return ((l - r) > NumErr<RealNum>::v);
}

template<typename Integer>
bool isOdd(Integer n) { return ((n & 1) == 1); }
template<typename Integer>
bool isEven(Integer n) { return ((n & 1) == 0); }

inline int8_t sign(int8_t i) { return ((i >> ((sizeof(i) * CHAR_BIT) - 2)) & -2) + 1; }
inline int16_t sign(int16_t i) { return ((i >> ((sizeof(i) * CHAR_BIT) - 2)) & -2) + 1; }
inline int32_t sign(int32_t i) { return ((i >> ((sizeof(i) * CHAR_BIT) - 2)) & -2) + 1; }
inline int64_t sign(int64_t i) { return ((i >> ((sizeof(i) * CHAR_BIT) - 2)) & -2) + 1; }
inline int32_t sign(float f) { return ((pCast<int32_t>(f) >> 30) & -2) + 1; }
inline int64_t sign(double f) { return ((pCast<int64_t>(f) >> 62) & -2) + 1; }

// one/both of them are 0: return false.
template<typename T>
bool sameSign(T l, T r) { return ((l > 0) == (r > 0)) || ((l < 0) == (r < 0)); }
//bool sameSign(T l, T r) { return (l * r) > 0; } // slow if `T` is `BigNumber`-like data structure.

// one/both of them is 0: return true.
template<typename T>
bool sameSign0(T l, T r) { return ((l >= 0) == (r >= 0)) || ((l <= 0) == (r <= 0)); }
//bool sameSign0(T l, T r) { return (l * r) >= 0; } // slow if `T` is `BigNumber`-like data structure.
inline bool sameSign0(short l, short r) { return ((l ^ r) >= 0) || (l == 0) || (r == 0); }
inline bool sameSign0(int l, int r) { return ((l ^ r) >= 0) || (l == 0) || (r == 0); }
inline bool sameSign0(long l, long r) { return ((l ^ r) >= 0) || (l == 0) || (r == 0); }
inline bool sameSign0(long long l, long long r) { return ((l ^ r) >= 0) || (l == 0) || (r == 0); }
inline bool sameSign0(float l, float r) { return ((pCast<int32_t>(l) ^ pCast<int32_t>(r)) >= 0) || (l == 0) || (r == 0); }
inline bool sameSign0(double l, double r) { return ((pCast<int64_t>(l) ^ pCast<int64_t>(r)) >= 0) || (l == 0) || (r == 0); }

template<typename T>
bool updateMin(T& minItem, const T& newItem) {
	if (newItem < minItem) {
		minItem = newItem;
		return true;
	}
	return false;
}
template<typename T, typename LessThan>
bool updateMin(T& minItem, const T& newItem, LessThan lessThan) {
	if (lessThan(newItem, minItem)) {
		minItem = newItem;
		return true;
	}
	return false;
}

template<typename T>
bool updateMax(T& maxItem, const T& newItem) {
	if (maxItem < newItem) {
		maxItem = newItem;
		return true;
	}
	return false;
}
template<typename T, typename LessThan>
bool updateMax(T& maxItem, const T& newItem, LessThan lessThan) {
	if (lessThan(maxItem, newItem)) {
		maxItem = newItem;
		return true;
	}
	return false;
}

template<typename T>
T bound(T num, T lb, T ub) {
	return (std::min)((std::max)(num, lb), ub);
}

template<typename T>
bool inRange(T num, T lb, T ub) { // check if `num \in [lb, ub)` is true.
	return (num >= lb) && (num < ub);
}

template<typename T>
bool between(T mid, T n0, T n1) { // check if `num \in [min(n0, n1), max(n0, n1)]` is true.
	return (mid - n0) * (mid - n1) < 0;
}
#pragma endregion Compare

#pragma region Statistics
template<typename T>
T average(T num1, T weight1, T num2, T weight2) {
	return ((num1 * weight1) + (num2 * weight2)) / (weight1 + weight2);
}

// normalize the list to [0, maxValue].
template<typename T, typename List>
void normalize(List& list, T maxValue, T minItem, T maxItem) {
	T range = maxItem - minItem;
	for (auto iter = list.begin(); iter != list.end(); ++iter) {
		((*iter -= minItem) *= maxValue) /= range;
	}
}
template<typename T, typename List>
void normalize(List& list, T maxValue) {
	auto itemRange = std::minmax_element(list.begin(), list.end());
	normalize(list, maxValue, *(itemRange.first), *(itemRange.second));
}
#pragma endregion Statistics

#pragma region Arithmetic
template<typename Integer = int>
struct GcdPositiveImpl { // for non-negative integers only.
	static Integer iterSub(Integer a, Integer b) {
		for (;;) {
			if (a > b) {
				a -= b;
			} else if (a < b) {
				b -= a;
			} else {
				return a;
			}
		}
	}

	static Integer recurDiv(Integer a, Integer b) { return (b == 0) ? a : recurDiv(b, a % b); }
	static Integer iterDiv(Integer a, Integer b) {
		for (;;) {
			if (b == 0) { return a; }
			a %= b;
			if (a == 0) { return b; }
			b %= a;
		}
	}

	static Integer recurStein(Integer a, Integer b) {
		if (a == 0) { return b; }
		if (b == 0) { return a; }
		if (isEven(a)) {
			if (isEven(b)) { return 2 * recurStein(a >> 1, b >> 1); }
			return recurStein(a >> 1, b);
		}
		if (isEven(b)) { return recurStein(a, b >> 1); }

		return (a > b) ? recurStein(a - b, b) : recurStein(b - a, a);
	}
	static Integer iterStein(Integer a, Integer b) {
		if (a == 0) { return b; }
		if (b == 0) { return a; }

		Integer k; // greatest power of 2 that divides both a and b.
		for (k = 0; isEven(a | b); ++k) {
			a >>= 1;
			b >>= 1;
		}

		while (isEven(a)) { a >>= 1; } // dividing `a` by 2 until `a` becomes odd.

		do { // from here on, `a` is always odd.
			while (isEven(b)) { b >>= 1; } // remove all factor of 2 in `b`.
			if (a > b) { std::swap(a, b); }
			b -= a;
		} while (b != 0);
		
		return a << k; // restore common factors of 2.
	}
};

template<typename Integer = int>
Integer gcd(Integer a, Integer b) { return GcdPositiveImpl<Integer>::iterDiv(a, b); }

// the result will always be a non-negative integer.
template<typename Integer = int>
Integer gcdp(Integer a, Integer b) { return GcdPositiveImpl<Integer>::iterStein(std::abs(a), std::abs(b)); }

template<typename Integer = int>
Integer lcm(Integer a, Integer b) { return (a / gcd(a, b)) * b; }

// the result of GCD will always be a non-negative integer.
template<typename Integer = int>
Integer lcmp(Integer a, Integer b) { return (a / gcdp(a, b)) * b; }

template<typename T>
constexpr T power2(T num) { return (num * num); }

class Log2Impl {
public:
	using Uint = uint32_t;

	static Uint v0(Uint n) {
		return sCast<Uint>(std::log2(n));
	}

	static Uint v1(Uint n);

	// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
	static Uint v2(Uint n) {
		#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
		static const char LogTable256[256] = {
			-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
			LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
			LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
		};
		#undef LT

		Uint r;
		Uint tt = (n >> 16);
		if (tt) {
			Uint t = (tt >> 8);
			r = (t ? (24 + LogTable256[t]) : (16 + LogTable256[tt]));
		} else {
			Uint t = (n >> 8);
			r = (t ? (8 + LogTable256[t]) : LogTable256[n]);
		}
		return r;
	}

	// http://guihaire.com/code/?p=414
	static Uint v3(Uint n) {
		union {
			float f;
			Uint n;
		} u;
		u.f = sCast<float>(n);
		return ((u.n >> 23) - 127);
	}

	// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogDeBruijn
	static Uint v4(Uint n) {
		static const Uint MultiplyDeBruijnBitPosition[32] = {
			0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
			8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
		};

		n |= n >> 1;
		n |= n >> 2;
		n |= n >> 4;
		n |= n >> 8;
		n |= n >> 16;

		return MultiplyDeBruijnBitPosition[sCast<Uint>(n * 0x07C4ACDDU) >> 27];
	}

	// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
	static Uint v5(Uint n) {
		const Uint b[] = { 0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000 };
		const Uint S[] = { 1, 2, 4, 8, 16 };

		Uint r = 0;
		for (int i = 4; i >= 0; --i) {
			if (n & b[i]) {
				n >>= S[i];
				r |= S[i];
			}
		}
		return r;
	}

	// http://stackoverflow.com/questions/9411823/fast-log2float-x-implementation-c/9411984#9411984
	static Uint v6(Uint n) {
		int r;
		std::frexp(n, &r);
		return sCast<Uint>(r - 1);
	}
};

// [Inaccurate]
inline Log2Impl::Uint log2(Log2Impl::Uint n) {
	#if _CC_MS_VC
	return Log2Impl::v1(n);
	#elif _CC_GNU_GCC
	return Log2Impl::v3(n);
	#else
	return Log2Impl::v2(n);
	#endif // _CC_MS_VC
}

// regard $0 == 2^{-\infty}$.
template<typename UnsignedIntegerType>
constexpr bool isPowerOf2(UnsignedIntegerType x) {
	return (x & (x - 1)) == 0;
}

template<typename IntegerType>
constexpr IntegerType leastSignificantBit(IntegerType i) {
	return i & (-i); // i & ((~i) + 1).
}

template<typename IntegerType>
IntegerType mostSignificantBit(IntegerType i) {
	i |= (i >> 1);
	i |= (i >> 2);
	i |= (i >> 4);
	if (sizeof(IntegerType) > 1) { i |= (i >> 8); }
	if (sizeof(IntegerType) > 2) { i |= (i >> 16); }
	if (sizeof(IntegerType) > 4) { i |= (i >> 32); }
	return i + 1;
}

template<typename IntegerType>
constexpr IntegerType lastBits(int n) { return (sCast<IntegerType>(1) << n) - 1; }

template<typename IntegerType>
bool isMultiplyOverflow(IntegerType l, IntegerType r) {
	if (r == 0) { return false; }
	IntegerType m = l * r;
	return (m / r) != l;
}

template<typename IntegerType, typename RealType>
IntegerType ceil(RealType x) { return sCast<IntegerType>(std::ceil(x)); }

template<typename IntegerType, typename RealType>
IntegerType floor(RealType x) { return sCast<IntegerType>(std::floor(x)); }

template<typename IntegerType>
IntegerType sqrtCeil(IntegerType i) { return ceil<IntegerType>(std::sqrt(i)); }

short sqrtFloor(short x);

int sqrtFloor(int x);

// https://stackoverflow.com/a/14100975
// Approximates atan2(y, x) normalized to the [0,4) range with a maximum error of 0.1620 degrees.
// atan(x) ~ `(c x + x^2 + x^3) / (1 + (c + 1) x + (c + 1) x^2 + x^3)` where `c = (1 + sqrt(17)) / 8` with a maximum error of 0.00811 degrees.
inline float arctan(float y, float x) {
	static constexpr uint32_t sign_mask = 0x80000000;
	static constexpr float b = 0.596227f;

	// Extract the sign bits
	uint32_t ux_s = sign_mask & (uint32_t&)x;
	uint32_t uy_s = sign_mask & (uint32_t&)y;

	// Determine the quadrant offset
	float q = (float)((~ux_s & uy_s) >> 29 | ux_s >> 30);

	// Calculate the arctangent in the first quadrant
	float bxy_a = ::fabs(b * x * y);
	float num = bxy_a + y * y;
	float atan_1q = num / (x * x + bxy_a + num + 1e-4f);

	// Translate it to the proper quadrant
	uint32_t uatan_2q = (ux_s ^ uy_s) | (uint32_t&)atan_1q;
	return q + (float&)uatan_2q;
}
#pragma endregion Arithmetic

}
}


#endif // CN_HUST_GOAL_COMMON_MATH_H
