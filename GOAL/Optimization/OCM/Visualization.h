////////////////////////////////
/// usage : 1.	visualize input and output by mathematica.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_O_C_M_VISUALIZATION_H
#define CN_HUST_GOAL_O_C_M_VISUALIZATION_H


#include <iostream>
#include <fstream>

#include "GOAL/Typedef.h"
#include "Problem.h"


namespace goal {
namespace OCM {

struct Visualizer {
	static void drawInput(const Str& path, const Input& input);
	static void drawOutput(const Str& path, const Input& input, const Output& output);

	struct Circle {
		Coord x, y;
	};
};

}
}


#endif // CN_HUST_GOAL_O_C_M_VISUALIZATION_H
