////////////////////////////////
/// usage : 1.	data structure for stable rational numbers.
/// 
/// note  : 1.	the denominator `div` is always positive.
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_RATIONAL_H
#define CN_HUST_GOAL_COMMON_RATIONAL_H


#include <numeric>
#include <functional>

#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Common/VarStr.h"
#include "GOAL/Math/Math.h"


namespace goal {
namespace math {

template<typename Integer = long long>
struct Rational {
	Rational() {}
	Rational(Integer numerator, Integer denominator) : num(numerator), div(denominator) { norm(); }
	Rational(double r) { fromReal(r); }

	size_t hash() const { return sCast<size_t>(num) * sCast<size_t>(div); }


	void normSign() {
		if (div < 0) { div = -div; num = -num; }
	}
	void cancel() {
		Integer c = gcd(std::abs(num), div);
		num /= c;
		div /= c;
	}
	void norm() {
		normSign();
		cancel();
	}

	void fromRealNaive(double r, Integer scale = (1 << 30)) {
		div = scale;
		num = sCast<Integer>(r * scale);
		cancel();
	}
	// https://stackoverflow.com/a/45314258
	// https://stackoverflow.com/a/37573546
	void fromRealBySjaak(double r, double err = NumErr<double>::v) {
		Integer intPart = (r > 0) ? sCast<Integer>(r) : -sCast<Integer>(1 - r);
		r -= intPart;

		double lb = r - err;
		if (lb < 0.0) {
			num = intPart;
			div = 1;
			return;
		}
		double ub = r + err;
		if (ub > 1.0) {
			num = intPart + 1;
			div = 1;
			return;
		}

		Integer b = 1;
		Integer d = sCast<Integer>(1 / ub);
		double lowNum = lb;
		double lowDiv = 1.0 - d * lb;
		double highNum = 1.0 - d * ub;
		double highDiv = ub;
		for (;;) {
			if (lowNum < lowDiv) { break; }
			Integer n = sCast<Integer>(lowNum / lowDiv);
			b += (n * d);
			lowNum -= (n * lowDiv);
			highDiv -= (n * highNum);
			if (highNum < highDiv) { break; }
			n = sCast<Integer>(highNum / highDiv);
			d += (n * b);
			lowDiv -= (n * lowNum);
			highNum -= (n * highDiv);
		}

		div = b + d;
		num = sCast<Integer>(r * div + 0.5) + (intPart * div);
	}
	void fromReal(double r) { fromRealBySjaak(r); }

	// num   r.num   num * (r.div / c)   r.num * (div / c)   num * lCoef + r.num * rCoef
	// --- + ----- = ----------------- + ----------------- = ---------------------------
	// div   r.div   div * (r.div / c)   r.div * (div / c)   div * lCoef = r.div * rCoef
	Rational& operator+=(const Rational& r) {
		Integer c = gcd(div, r.div);
		Integer lCoef = r.div / c;
		Integer rCoef = div / c;
		div *= lCoef;
		num *= lCoef;
		num += (r.num * rCoef);
		norm();
		return *this;
	}
	// num   r.num   num * (r.div / c)   r.num * (div / c)   num * lCoef - r.num * rCoef
	// --- - ----- = ----------------- - ----------------- = ---------------------------
	// div   r.div   div * (r.div / c)   r.div * (div / c)   div * lCoef = r.div * rCoef
	//Rational& operator-=(const Rational& r) { return operaotr+=(Rational(-r.num, r.div)); }
	Rational& operator-=(const Rational& r) {
		Integer c = gcd(div, r.div);
		Integer lCoef = r.div / c;
		Integer rCoef = div / c;
		div *= lCoef;
		num *= lCoef;
		num -= (r.num * rCoef);
		norm();
		return *this;
	}
	// num   r.num   num / nd   r.num / dn
	// --- * ----- = -------- * ----------
	// div   r.div   div / dn   r.div / nd
	Rational& operator*=(const Rational& r) {
		Integer nd = gcd(std::abs(num), r.div);
		Integer dn = gcd(div, std::abs(r.num));
		(num /= nd) *= (r.num / dn);
		(div /= dn) *= (r.div / nd);
		return *this;
	}
	// num   r.num   num / nn   r.div / dd
	// --- / ----- = -------- * ----------
	// div   r.div   div / dd   r.num / nn
	//Rational& operator/=(const Rational& r) { return operaotr*=(Rational(r.div, r.num)); }
	Rational& operator/=(const Rational& r) {
		Integer nn = gcd(std::abs(num), std::abs(r.num));
		Integer dd = gcd(div, r.div);
		(num /= nn) *= (r.div / dd);
		(div /= dd) *= (r.num / nn);
		normSign();
		return *this;
	}

	friend Rational operator+(const Rational& l, const Rational& r) { return Rational(l) += r; }
	friend Rational operator-(const Rational& l, const Rational& r) { return Rational(l) -= r; }
	friend Rational operator*(const Rational& l, const Rational& r) { return Rational(l) *= r; }
	friend Rational operator/(const Rational& l, const Rational& r) { return Rational(l) /= r; }

	Rational& operator++() { num += div; return *this; }
	Rational& operator--() { num -= div; return *this; }
	Rational operator++(int) {
		Rational r(*this);
		++(*this);
		return r;
	}
	Rational operator--(int) {
		Rational r(*this);
		--(*this);
		return r;
	}

	Rational operator-() const { return Rational(-num, div); }

	friend bool operator==(const Rational& l, const Rational& r) { return (l.num == r.num) && (l.div == r.div); }
	friend bool operator!=(const Rational& l, const Rational& r) { return !(l == r); }
	friend bool operator<(const Rational& l, const Rational& r) {
		Integer c = gcd(l.div, r.div);
		Integer lCoef = r.div / c;
		Integer rCoef = l.div / c;
		return (l.num * lCoef) < (r.num * rCoef);
	}
	friend bool operator>(const Rational& l, const Rational& r) {
		Integer c = gcd(l.div, r.div);
		Integer lCoef = r.div / c;
		Integer rCoef = l.div / c;
		return (l.num * lCoef) > (r.num * rCoef);
	}
	friend bool operator<=(const Rational& l, const Rational& r) { return !(l > r); }
	friend bool operator>=(const Rational& l, const Rational& r) { return !(l < r); }

	//template<typename OutputStream>
	//friend OutputStream& operator<<(OutputStream& os, const Rational& r) { return os << r.num << '/' << r.div; }
	//template<typename InputStream>
	//friend InputStream& operator>>(InputStream& is, Rational& r) { char c; return is >> num >> c >> div; }
	template<typename OutputStream>
	friend OutputStream& operator<<(OutputStream& os, const Rational& r) { return os << sCast<double>(r); }
	template<typename InputStream>
	friend InputStream& operator>>(InputStream& is, Rational& r) {
		double d;
		is >> d;
		r = Rational(d);
		return is;
	}

	explicit operator double() const { return sCast<double>(num) / sCast<double>(div); }
	operator Str() const { return toStr(num) + '/' + toStr(div); }

protected:
	Integer num;
	Integer div;
};

}
}

template<typename Integer>
struct std::hash<goal::math::Rational<Integer>> {
	size_t operator()(goal::math::Rational<Integer> const& r) const { return r.hash(); }
};


#endif // CN_HUST_GOAL_COMMON_RATIONAL_H
