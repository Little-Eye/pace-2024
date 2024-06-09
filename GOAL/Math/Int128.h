////////////////////////////////
/// usage : 1.	data structure for 128-bit signed integer.
/// 
/// note  : 1.	i[0] is unsigned. i[1] is signed.
////////////////////////////////

#ifndef CN_HUST_GOAL_MATH_INT128_H
#define CN_HUST_GOAL_MATH_INT128_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Math/Math.h"

#include <cstdint>


namespace goal {
namespace math {

struct Uint128 {
	using Uint32 = uint32_t;
	using Uint64 = uint64_t;


	static constexpr int QuadBytes = sizeof(Uint32);
	static constexpr int QuadBits = QuadBytes * 8;
	static constexpr int HalfBytes = sizeof(Uint64);
	static constexpr int HalfBits = HalfBytes * 8;
	static constexpr int FullBytes = HalfBytes * 2;
	static constexpr int FullBits = FullBytes * 8;
	static constexpr int SignBitIndex = HalfBits - 1; // zero-based.
	static constexpr Uint64 LowMask = sCast<Uint32>(-1);
	static constexpr Uint64 HighMask = (LowMask << QuadBits);


	Uint128() {}
	Uint128(Uint64 l, Uint64 h = 0) { i[0] = l; i[1] = h; }
	Uint128(Arr<Uint64, 2> n) { i = n; }


	Uint128& operator&=(const Uint128& r) {
		i[1] &= r.i[1];
		i[0] &= r.i[0];
		return *this;
	}
	Uint128& operator|=(const Uint128& r) {
		i[1] |= r.i[1];
		i[0] |= r.i[0];
		return *this;
	}
	Uint128& operator^=(const Uint128& r) {
		i[1] ^= r.i[1];
		i[0] ^= r.i[0];
		return *this;
	}
	Uint128& operator>>=(int r) {
		Uint64 mask = lastBits<Uint64>(r);
		i[0] >>= r;
		i[0] |= ((i[1] & mask) << (HalfBits - r));
		i[1] >>= r;
		return *this;
	}
	Uint128& operator<<=(int r) {
		Uint64 mask = lastBits<Uint64>(r);
		i[1] <<= r;
		i[1] |= ((i[0] >> (HalfBits - r)) & mask);
		i[0] <<= r;
		return *this;
	}
	Uint128& shiftRight1() {
		i[0] >>= 1;
		i[0] |= (i[1] << (HalfBits - 1));
		i[1] >>= 1;
		return *this;
	}
	Uint128& shiftLeft1() {
		i[1] <<= 1;
		i[1] |= (i[0] >> (HalfBits - 1));
		i[0] <<= 1;
		return *this;
	}
	static Uint128 shiftRight1(Uint128 n) { return n.shiftRight1(); }
	static Uint128 shiftLeft1(Uint128 n) { return n.shiftLeft1(); }

	Uint128& operator+=(const Uint128& r) {
		i[1] += r.i[1] + ((i[0] + r.i[0]) < i[0]);
		i[0] += r.i[0];
		return *this;
	}
	Uint128& operator-=(const Uint128& r) {
		i[1] -= (r.i[1] + ((i[0] - r.i[0]) > i[0]));
		i[0] -= r.i[0];
		return *this;
	}
	Uint128& operator*=(const Uint128& r) {
		if (r == 0u) { i[0] = i[1] = 0; return *this; }
		if (r == Uint128(-1, -1)) { return *this = -*this; }

		Uint128 e[FullBits];
		e[0] = *this;
		int i = 0;
		Uint128 k = 1;
		for (; k.shiftLeft1() < r; ++i) { e[i + 1] = shiftLeft1(e[i]); }

		*this = e[i];
		Uint128 rr = r;
		while ((rr -= k.shiftRight1()) > 0u) {
			for (k = 1, i = 0; k < rr; k.shiftLeft1(), ++i) {}
			*this += e[i];
		}

		return *this;
	}
	Uint128& operator*=(const Uint128& r) {
		if (r == 0u) { i[0] = i[1] = 0; return *this; }

		Uint128 e[FullBits];
		e[0] = *this;
		int i = 0;
		Uint128 k = 1;
		for (; k < r; k.shiftLeft1()) { e[++i] = shiftLeft1(e[i]); }

		*this = 0;
		Uint128 r0 = r - k;
		Uint128 r1 = shiftLeft1(k) - r;
		if (r0 <= r1) {
			*this += e[i];
			for (k = 1, i = 0; k < r0; k.shiftLeft1(), ++i) {}
		} else {
			*this += shiftLeft1(e[i]);
			for (k = 1, i = 0; k < r0; k.shiftLeft1(), ++i) {}
			*this -= e[i];
		}

		return *this;
	}
	Uint128& operator*=(const Uint128& r) {
		//        x3 x2 x1 x0
		//   x    y3 y2 y1 y0
		//	 ---------------------
		//      | 30 20 10 00
		// over | 21 11 01
		// flow | 12 02
		//      | 33
		const Uint32* pl = sCast<const Uint32*>(i.data());
		const Uint32* pr = sCast<const Uint32*>(r.i.data());
		Uint64 x[QuadBytes] = { pl[0], pl[1], pl[2], pl[3] };
		Uint64 y[QuadBytes] = { pr[0], pr[1], pr[2], pr[3] };

		i[0] = x[0] * y[0];
		i[1] = x[1] * y[0] + (i[0] & HighMask);
		i[2] = x[2] * y[0];
		i[3] = x[3] * y[0];

		return *this;
	}
	Uint128& operator/=(const Uint128& r) {
		//num /= r.num;
		return *this;
	}

	friend Uint128 operator+(const Uint128& l, const Uint128& r) { return Uint128(l) += r; }
	friend Uint128 operator-(const Uint128& l, const Uint128& r) { return Uint128(l) -= r; }
	friend Uint128 operator*(const Uint128& l, const Uint128& r) { return Uint128(l) *= r; }
	friend Uint128 operator/(const Uint128& l, const Uint128& r) { return Uint128(l) /= r; }

	Uint128& operator++() { if (++i[0] == 0) { ++i[1]; } return *this; }
	Uint128& operator--() { if (--i[0] == -1) { --i[1]; }; return *this; }
	//Uint128 operator++(int) { return Uint128(num++); }
	//Uint128 operator--(int) { return Uint128(num--); }
	static Uint128 inc(Uint128 n) { return ++n; }
	static Uint128 dec(Uint128 n) { return --n; }

	Uint128 operator~() const { return Uint128(~i[0], ~i[1]); }
	Uint128 operator-() const { return inc(~*this); }

	friend bool operator==(const Uint128& l, const Uint128& r) { return l.i == r.i; }
	friend bool operator!=(const Uint128& l, const Uint128& r) { return l.i != r.i; }
	friend bool operator<(const Uint128& l, const Uint128& r) { return (l.i[1] < r.i[1]) || ((l.i[1] == r.i[1]) && (l.i[0] < r.i[0])); }
	friend bool operator>(const Uint128& l, const Uint128& r) { return (l.i[1] > r.i[1]) || ((l.i[1] == r.i[1]) && (l.i[0] > r.i[0])); }
	friend bool operator<=(const Uint128& l, const Uint128& r) { return (l.i[1] < r.i[1]) || ((l.i[1] == r.i[1]) && (l.i[0] <= r.i[0])); }
	friend bool operator>=(const Uint128& l, const Uint128& r) { return (l.i[1] > r.i[1]) || ((l.i[1] == r.i[1]) && (l.i[0] >= r.i[0])); }

	friend bool operator==(const Uint128& l, Uint64 r) { return (l.i[0] == r) && (l.i[1] == 0); }
	friend bool operator!=(const Uint128& l, Uint64 r) { return (l.i[0] != r) || (l.i[1] != 0); }
	friend bool operator<(const Uint128& l, Uint64 r) { return (l.i[1] == 0) && (l.i[0] < r); }
	friend bool operator>(const Uint128& l, Uint64 r) { return (l.i[1] > 0) || (l.i[0] > r); }
	friend bool operator<=(const Uint128& l, Uint64 r) { return (l.i[1] == 0) && (l.i[0] <= r); }
	friend bool operator>=(const Uint128& l, Uint64 r) { return (l.i[1] > 0) || (l.i[0] >= r); }

	friend bool operator==(const Uint128& l, Uint32 r) { return (l.i[0] == r) && (l.i[1] == 0); }
	friend bool operator!=(const Uint128& l, Uint32 r) { return (l.i[0] != r) || (l.i[1] != 0); }
	friend bool operator<(const Uint128& l, Uint32 r) { return (l.i[1] == 0) && (l.i[0] < r); }
	friend bool operator>(const Uint128& l, Uint32 r) { return (l.i[1] > 0) || (l.i[0] > r); }
	friend bool operator<=(const Uint128& l, Uint32 r) { return (l.i[1] == 0) && (l.i[0] <= r); }
	friend bool operator>=(const Uint128& l, Uint32 r) { return (l.i[1] > 0) || (l.i[0] >= r); }

	//template<typename OutputStream>
	//friend OutputStream& operator<<(OutputStream& os, const Uint128& r) { return os << r.num; }
	//template<typename InputStream>
	//friend InputStream& operator>>(InputStream& is, const Uint128& r) { return is >> r.num; }

	//explicit operator Rep() const { return num; }
	//operator Str() const { return toStr(num); }


protected:
	Arr<Uint64, 2> i;
};

struct Int128 {
	using Int64 = int64_t;
	using Uint64 = uint64_t;


	static constexpr int Int32ByteNum = sizeof(uint32_t);
	static constexpr int SignBitIndex64 = 64 - 1; // zero-based.
	static constexpr int High32Index64 = 32; // zero-based.
	static constexpr Int64 SignMask64 = (1ll << SignBitIndex64);
	static constexpr Int64 Low32Mask64 = 0xFFFFFFFFll;
	static constexpr Int64 High32Mask64 = (Low32Mask64 << High32Index64);


	Int128() {}
	Int128(Int64 n) { i[0] = n; i[1] = (n & SignMask64) ? -1 : 0; }
	Int128(Int64 h, Uint64 l) { i[0] = l; i[1] = h; }
	Int128(Arr<Int64, 2> n) { i = n; }

	Int128& operator+=(Int128 r) {
		i[1] += ((i[0] & r.i[0]) >> SignBitIndex64);
		i[1] += r.i[1];
		i[0] += r.i[0];
		return *this;
	}
	Int128& operator-=(Int128 r) {
		i[1] -= ((i[0] - r.i[0]) >> SignBitIndex64);
		i[1] -= r.i[1];
		i[0] -= r.i[0];
		return *this;
	}
	Int128& operator*=(Int128 r) {
		//         a  b  c  d
		//   x     e  f  g  h
		//	 -----------------
		//      | ah bh ch dh
		// over | bg cg dg
		// flow | cf df
		//      | de
		Uint64 x[Int32ByteNum];
		Uint64 y[Int32ByteNum];
		uint32_t* pl = sCast<uint32_t*>(i);
		uint32_t* pr = sCast<uint32_t*>(r.i);
		for (int k = 0; k < Int32ByteNum; ++k) { x[k] = pl[k]; y[k] = pr[k]; }


		Uint64 acc, ac2, carry;
		Uint64 a, b, c, d, e, f, g, h;
		d = i[0] & Low32Mask64;
		c = (i[0] & High32Mask64) >> High32Index64;
		b = i[1] & Low32Mask64;
		a = (i[1] & High32Mask64) >> High32Index64;

		h = r.i[0] & Low32Mask64;
		g = (r.i[0] & High32Mask64) >> High32Index64;
		f = r.i[1] & Low32Mask64;
		e = (r.i[1] & High32Mask64) >> High32Index64;

		acc = d * h;
		i[0] = acc & Low32Mask64;
		acc >>= High32Index64;
		carry = 0;
		ac2 = acc + c * h; if (ac2 < acc) { carry++; }
		acc = ac2 + d * g; if (acc < ac2) { carry++; }
		i[0] |= (acc << High32Index64);
		ac2 = (acc >> High32Index64) | (carry << High32Index64);
		carry = 0;

		acc = ac2 + b * h; if (acc < ac2) { carry++; }
		ac2 = acc + c * g; if (ac2 < acc) { carry++; }
		acc = ac2 + d * f; if (acc < ac2) { carry++; }
		i[1] = acc & Low32Mask64;
		ac2 = (acc >> High32Index64) | (carry << High32Index64);

		acc = ac2 + a * h;
		ac2 = acc + b * g;
		acc = ac2 + c * f;
		ac2 = acc + d * e;
		i[1] |= (ac2 << High32Index64);
		return *this;
	}
	Int128& operator/=(Int128 r) {
		//num /= r.num;
		return *this;
	}

	friend Int128 operator+(const Int128& l, const Int128& r) { return Int128(l) += r; }
	friend Int128 operator-(const Int128& l, const Int128& r) { return Int128(l) -= r; }
	friend Int128 operator*(const Int128& l, const Int128& r) { return Int128(l) *= r; }
	friend Int128 operator/(const Int128& l, const Int128& r) { return Int128(l) /= r; }

	//Int128& operator++() { ++num; return *this; }
	//Int128& operator--() { --num; return *this; }
	//Int128 operator++(int) { return Int128(num++); }
	//Int128 operator--(int) { return Int128(num--); }

	Int128 operator-() const {
		// todo




		return Int128(-i[1], i[0]);
	}

	//friend bool operator==(Int128 l, Int128 r) { return math::weakEqual(l.num, r.num); }
	//friend bool operator!=(Int128 l, Int128 r) { return !math::weakEqual(l.num, r.num); }
	//friend bool operator<(Int128 l, Int128 r) { return math::strictLess(l.num, r.num); }
	//friend bool operator>(Int128 l, Int128 r) { return math::strictGreater(l.num, r.num); }
	//friend bool operator<=(Int128 l, Int128 r) { return math::weakLess(l.num, r.num); }
	//friend bool operator>=(Int128 l, Int128 r) { return math::weakGreater(l.num, r.num); }

	//template<typename OutputStream>
	//friend OutputStream& operator<<(OutputStream& os, Int128 r) { return os << r.num; }
	//template<typename InputStream>
	//friend InputStream& operator>>(InputStream& is, Int128 r) { return is >> r.num; }

	//explicit operator Rep() const { return num; }
	//operator Str() const { return toStr(num); }

protected:
	Arr<Int64, 2> i;
};

}
}


#endif // CN_HUST_GOAL_MATH_INT128_H
