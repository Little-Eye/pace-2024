////////////////////////////////
/// usage : 1.	data structures and algorithms for computational geometry.
/// 
/// note  : 1.	it is better to use `math::SafeReal` or `math::Rational` as template parameter `Number` to overcome precision issues.
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_GEOMETRY_H
#define CN_HUST_GOAL_COMMON_GEOMETRY_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/Graph/Graph.h"
#include "GOAL/Math/Math.h"


namespace goal {
namespace graph {

template<typename Number = Real>
using Point2D = Coord2D<Number>;

template<typename Number = Real>
struct Line2D : public Arr<Point2D<Number>, 2> {
	Line2D(const Point2D<Number>& p0, const Point2D<Number>& p1) : Arr<Point2D<Number>, 2>({ p0, p1 }) {}
};

template<typename Number = Real>
struct Ray2D : public Line2D<Number> { using Line2D<Number>::Line2D; }; // heading from `Ray2D[0]` to `Ray2D[1]`.

template<typename Number = Real>
struct Segment2D : public Line2D<Number> { using Line2D<Number>::Line2D; };

template<typename Number = Real>
struct Angle2D : public Arr<Point2D<Number>, 3> { // heading from `Angle2D[0]` to `Angle2D[1]` to `Angle2D[2]`.
	Angle2D(const Point2D<Number>& p0, const Point2D<Number>& p1, const Point2D<Number>& p2) : Arr<Point2D<Number>, 3>({ p0, p1, p2 }) {}
};


inline double polarRay(double x0, double y0, double x1, double y1) { // the angle between segment `s` and ray (0, 0)->(0, 1).
	return std::atan2(x1 - x0, y1 - y0); // https://en.cppreference.com/w/cpp/numeric/math/atan2
}

template<typename Number = Real>
double polarRay(const Point2D<Number>& p0, const Point2D<Number>& p1) {
	return polarRay(sCast<double>(p0[0]), sCast<double>(p0[1]), sCast<double>(p1[0]), sCast<double>(p1[1]));
}

template<typename Number = Real>
void lineEquation(Number& xCoef, Number& yCoef, Number& constant, Number x0, Number y0, Number x1, Number y1) {
	xCoef = y1 - y0; // A.
	yCoef = x0 - x1; // B.
	constant = x1 * y0 - x0 * y1; // C.
}

template<typename Number = Real>
void lineEquation(Number& xCoef, Number& yCoef, Number& constant, const Line2D<Number>& line) {
	lineEquation(xCoef, yCoef, constant, line[0][0], line[0][1], line[1][0], line[1][1]);
}


template<typename Number = Real>
bool intersectLL(Point2D<Number>& intersection, const Line2D<Number>& line0, const Line2D<Number>& line1) {
	Number a0, b0, c0, a1, b1, c1;
	lineEquation(a0, b0, c0, line0);
	lineEquation(a1, b1, c1, line1);
	Number denominator = a0 * b1 - a1 * b0;
	if (denominator == 0) { return false; } // FIX[szx][0]: what if they overlap???!!!
	intersection[0] = (c1 * b0 - c0 * b1) / denominator;
	intersection[1] = (c0 * a1 - c1 * a0) / denominator;
	return true;
}
template<typename Number = Real>
bool intersect(Point2D<Number>& intersection, const Line2D<Number>& line0, const Line2D<Number>& line1) {
	return intersectLL(intersection, line0, line1);
}

template<typename Number = Real>
bool intersectRL(Point2D<Number>& intersection, const Ray2D<Number>& ray, const Line2D<Number>& line) {
	if (!intersectLL(intersection, ray, line)) { return false; }
	Number dx = ray[1][0] - ray[0][0];
	if (dx != 0) {
		return math::sameSign0(intersection[0] - ray[0][0], dx);
	} else { // vertical ray.
		return math::sameSign0(intersection[1] - ray[0][1], ray[1][1] - ray[0][1]);
	}
}
template<typename Number = Real>
bool intersect(Point2D<Number>& intersection, const Ray2D<Number>& ray, const Line2D<Number>& line) {
	return intersectRL(intersection, ray, line);
}

template<typename Number = Real>
bool intersectRS(Point2D<Number>& intersection, const Ray2D<Number>& ray, const Segment2D<Number>& seg) {
	if (!intersectRL(intersection, ray, seg)) { return false; }
	Number dx = seg[1][0] - seg[0][0];
	if (dx != 0) {
		return math::sameSign0(intersection[0] - seg[0][0], seg[1][0] - intersection[0]);
	} else { // vertical segment.
		return math::sameSign0(intersection[1] - seg[0][1], seg[1][1] - intersection[1]);
	}
}
template<typename Number = Real>
bool intersect(Point2D<Number>& intersection, const Ray2D<Number>& ray, const Segment2D<Number>& seg) {
	return intersectRS(intersection, ray, seg);
}

template<typename Number = Real>
bool intersect(Point2D<Number>& intersection, const Segment2D<Number>& seg0, const Segment2D<Number>& seg1) {
	// FIX[szx][5]: intersect
}


template<typename Number = Real>
bool isOnLeft(Number x0, Number y0, Number x1, Number y1, Number x, Number y) {
	return (x1 - x0) * (y - y0) > (y1 - y0) * (x - x0);
}

template<typename Number = Real>
bool isOnLeft(const Ray2D<Number>& ray, const Point2D<Number>& point) {
	return isOnLeft<Number>(ray[0][0], ray[0][1], ray[1][0], ray[1][1], point[0], point[1]);
}

template<typename Number = Real>
bool isOnRight(Number x0, Number y0, Number x1, Number y1, Number x, Number y) {
	return (x1 - x0) * (y - y0) < (y1 - y0) * (x - x0);
}

template<typename Number = Real>
bool isOnRight(const Ray2D<Number>& ray, const Point2D<Number>& point) {
	return isOnRight<Number>(ray[0][0], ray[0][1], ray[1][0], ray[1][1], point[0], point[1]);
}

template<typename Number = Real>
bool isOnLine(const Ray2D<Number>& ray, const Point2D<Number>& point) {
	return (ray[1][0] - ray[0][0]) * (point[1] - ray[0][1]) == (ray[1][1] - ray[0][1]) * (point[0] - ray[0][0]);
}

template<typename Number = Real>
bool isOnSegment(const Ray2D<Number>& ray, const Point2D<Number>& point) {
	// FIX[szx][5]: isOnSegment
}

template<typename Number = Real>
bool isOnRay(const Ray2D<Number>& ray, const Point2D<Number>& point) {
	// FIX[szx][5]: isOnRay
}

template<typename Number = Real>
bool isTurnLeft(const Angle2D<Number>& angle) {
	return isOnLeft<Number>(angle[0][0], angle[0][1], angle[1][0], angle[1][1], angle[2][0], angle[2][1]);
}

template<typename Number = Real>
bool isTurnRight(const Angle2D<Number>& angle) {
	return isOnRight<Number>(angle[0][0], angle[0][1], angle[1][0], angle[1][1], angle[2][0], angle[2][1]);
}

template<typename Number = Real>
bool isGoForward(const Angle2D<Number>& angle) {
	return isOnRay<Number>(angle[0][0], angle[0][1], angle[1][0], angle[1][1], angle[2][0], angle[2][1]);
}

template<typename Number = Real>
bool isGoBackward(const Angle2D<Number>& angle) {
	return isOnRay<Number>(angle[0][0], angle[0][1], angle[1][0], angle[1][1], angle[2][0], angle[2][1]);
}

template<typename Number = Real>
Number area(const Vec<Point2D<Number>>& boundaryPoints) { // Shoelace Theorem: https://zhuanlan.zhihu.com/p/110025234
	ID num = sCast<ID>(boundaryPoints.size());
	Number a = (boundaryPoints[num - 1][0] * boundaryPoints[0][1]) - (boundaryPoints[0][0] * boundaryPoints[num - 1][1]);
	for (ID i = 1; i < num; ++i) {
		a += (boundaryPoints[i - 1][0] * boundaryPoints[i][1]);
		a -= (boundaryPoints[i][0] * boundaryPoints[i - 1][1]);
	}
	if (a < 0) { a = -a; }
	return a /= 2;
}

}
}


#endif // CN_HUST_GOAL_COMMON_GEOMETRY_H
