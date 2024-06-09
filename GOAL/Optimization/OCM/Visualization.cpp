#include "Visualization.h"

#include <iostream>
#include <sstream>
#include <iomanip>

#include "GOAL/Common/Log.h"
#include "GOAL/Common/Math.h"
#include "GOAL/System/SvgDrawer.h"


using namespace std;


namespace goal {
namespace OCM {

void Visualizer::drawInput(const Str& path, const Input& input) {
	// Box b;
	// SvgDrawer<Coord> sd(drawGraph(b, input));
	// ofstream ofs(path + ".html");
	// ofs << sd.toStr(b.minX - 10, b.minY - 10, b.maxX - b.minX + 20, b.maxY - b.minY + 20, 1920, 1080, "", "",
	// 	"<style>.b, .h { stroke:black; fill:transparent; }</style>") << endl;
	/*SvgDrawer<Coord> sd;
	Coord x, y, r = 1;
	Coord centerGap = 3;
	Vec<Circle> circle1(input.nodeNum1), circle2(input.nodeNum2);
	Coord viewBoxWidth = 1050, viewBoxHeight = 0;

	x = 1, y = 1;
	for (Index index = 0; index < input.nodeNum1; ++index) {
		NodeId nodeId = index;
		circle1[nodeId].x = x;
		circle1[nodeId].y = y;
		y += centerGap;
	}
	viewBoxHeight = max(viewBoxWidth, y);

	x = 1000, y = 1;
	for (Index index = 0; index < input.nodeNum2; ++index) {
		NodeId nodeId = index;
		circle2[nodeId].x = x;
		circle2[nodeId].y = y;
		y += centerGap;
	}
	viewBoxHeight = max(viewBoxWidth, y);

	for (NodeId nodeId = 0; nodeId < input.nodeNum1; ++nodeId) {
		sd.circle(circle1[nodeId].x, circle1[nodeId].y, r);
	}
	for (NodeId nodeId = 0; nodeId < input.nodeNum2; ++nodeId) {
		sd.circle(circle2[nodeId].x, circle2[nodeId].y, r);
	}
	for (NodeId nodeId2 = 0; nodeId2 < input.nodeNum2; ++nodeId2) {
		for (auto nodeId1 : input.inNodeIds[nodeId2]) {
			sd.line(circle1[nodeId1].x, circle1[nodeId1].y, circle2[nodeId2].x, circle2[nodeId2].y);
		}
	}

	ofstream ofs(path + ".html");
	ofs << sd.toStr(0, 0, viewBoxWidth, viewBoxHeight);*/
}

void Visualizer::drawOutput(const Str& path, const Input& input, const Output& output) {
	// if (output.convexPolygons.empty()) { return; }
	/*SvgDrawer<Coord> sd;
	Coord x, y, r = 1;
	Coord centerGap = 3;
	Vec<Circle> circle1(input.nodeNum1), circle2(input.nodeNum2);
	Coord viewBoxWidth = 1050, viewBoxHeight = 0;

	x = 1, y = 1;
	for (Index index = 0; index < input.nodeNum1; ++index) {
		NodeId nodeId = index;
		circle1[nodeId].x = x;
		circle1[nodeId].y = y;
		y += centerGap;
	}
	viewBoxHeight = max(viewBoxWidth, y);

	x = 1000, y = 1;
	for (Index index = 0; index < input.nodeNum2; ++index) {
		NodeId nodeId = output.order[index];
		circle2[nodeId].x = x;
		circle2[nodeId].y = y;
		y += centerGap;
	}
	viewBoxHeight = max(viewBoxWidth, y);

	for (NodeId nodeId = 0; nodeId < input.nodeNum1; ++nodeId) {
		sd.circle(circle1[nodeId].x, circle1[nodeId].y, r);
	}
	for (NodeId nodeId = 0; nodeId < input.nodeNum2; ++nodeId) {
		sd.circle(circle2[nodeId].x, circle2[nodeId].y, r);
	}
	for (NodeId nodeId2 = 0; nodeId2 < input.nodeNum2; ++nodeId2) {
		for (auto nodeId1 : input.inNodeIds[nodeId2]) {
			sd.line(circle1[nodeId1].x, circle1[nodeId1].y, circle2[nodeId2].x, circle2[nodeId2].y);
		}
	}

	ofstream ofs(path + ".html");
	ofs << sd.toStr(0, 0, viewBoxWidth, viewBoxHeight);*/
}

}
}
