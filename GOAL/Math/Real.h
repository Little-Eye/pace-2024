////////////////////////////////
/// usage : 1.	data structure for safe real numbers.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_MATH_REAL_H
#define CN_HUST_GOAL_MATH_REAL_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Common/VarStr.h"
#include "GOAL/Math/Math.h"


namespace goal {
namespace math {

template<typename Rep = double>
struct SafeReal {
	SafeReal() {}
	SafeReal(Rep realNum) : num(realNum) {}

	SafeReal& operator+=(SafeReal r) { num += r.num; return *this; }
	SafeReal& operator-=(SafeReal r) { num -= r.num; return *this; }
	SafeReal& operator*=(SafeReal r) { num *= r.num; return *this; }
	SafeReal& operator/=(SafeReal r) { num /= r.num; return *this; }

	friend SafeReal operator+(SafeReal l, SafeReal r) { return l.num + r.num; }
	friend SafeReal operator-(SafeReal l, SafeReal r) { return l.num - r.num; }
	friend SafeReal operator*(SafeReal l, SafeReal r) { return l.num * r.num; }
	friend SafeReal operator/(SafeReal l, SafeReal r) { return l.num / r.num; }

	SafeReal& operator++() { ++num; return *this; }
	SafeReal& operator--() { --num; return *this; }
	SafeReal operator++(int) { return SafeReal(num++); }
	SafeReal operator--(int) { return SafeReal(num--); }

	SafeReal operator-() const { return SafeReal(-num); }

	friend bool operator==(SafeReal l, SafeReal r) { return math::weakEqual(l.num, r.num); }
	friend bool operator!=(SafeReal l, SafeReal r) { return !math::weakEqual(l.num, r.num); }
	friend bool operator<(SafeReal l, SafeReal r) { return math::strictLess(l.num, r.num); }
	friend bool operator>(SafeReal l, SafeReal r) { return math::strictGreater(l.num, r.num); }
	friend bool operator<=(SafeReal l, SafeReal r) { return math::weakLess(l.num, r.num); }
	friend bool operator>=(SafeReal l, SafeReal r) { return math::weakGreater(l.num, r.num); }

	template<typename OutputStream>
	friend OutputStream& operator<<(OutputStream& os, SafeReal r) { return os << r.num; }
	template<typename InputStream>
	friend InputStream& operator>>(InputStream& is, SafeReal r) { return is >> r.num; }

	explicit operator Rep() const { return num; }
	operator Str() const { return toStr(num); }

protected:
	Rep num;
};

}
}


#endif // CN_HUST_GOAL_MATH_REAL_H
