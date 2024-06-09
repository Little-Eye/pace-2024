////////////////////////////////
/// usage : 1.	basic utilities for CG_DESCENT solver.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_THIRD_PARTY_CG_DESCENT_H
#define CN_HUST_GOAL_THIRD_PARTY_CG_DESCENT_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"

//#if _PLUGIN_CG_DESCENT

#include "CgDescent/cg_descent.h"


namespace goal {

class MpSolverCgDescent {
public:
	#pragma region Constant
	// status for the most recent optimization.
	enum StatusFlag {
		Optimal = 1,                     // CG_ERROR_TOLERANCE_SATISFIED
		SubOptimal = 1 << 1,             // CG_ITERATIONS_EXCEED_MAXITS || CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS
		FeasibleMask = Optimal | SubOptimal,

		InsolubleModel = 1 << 2,         // CG_SLOPE_ALWAYS_NEGATIVE
		IllConditioned = 1 << 3,         // CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION || CG_WOLFE_CONDITIONS_NOT_SATISFIED || CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES || CG_NO_COST_OR_GRADIENT_IMPROVEMENT
		ModelErrorMask = InsolubleModel | IllConditioned,

		OutOfMemory = 1 << 8,            // OUT_OF_MEMORY
		Unknown = 1 << 9,                // any other status code or error code.
		SolverErrorMask = OutOfMemory | Unknown,
	};
	using StatusFlags = int;
	#pragma endregion Constant

	#pragma region Type
	using Model = CGdata*;
	using Decision = CGFLOAT;
	using Index = CGINT;
	#pragma endregion Type

	#pragma region Constructor
	MpSolverCgDescent() : model(cg_setup()) {}

	~MpSolverCgDescent() { cg_terminate(&model); }
	#pragma endregion Constructor

	#pragma region Method
	StatusFlags optimize() {
		model->n = sCast<Index>(x.size());
		model->x = x.data();

		resetStatus();
		cg_descent(model);
		return getStatus();
	}

	StatusFlags optimize(CgEval value, CgEval grad) {
		setEvaluation(value, grad);
		return optimize();
	}
	StatusFlags optimize(CgEval value, CgEval grad, CgEval2 valgrad) {
		setEvaluation(value, grad, valgrad);
		return optimize();
	}

	StatusFlags getStatus() {
		status &= StatusFlag::SolverErrorMask;

		switch (model->Stat->status) {
		case CG_ERROR_TOLERANCE_SATISFIED:
			status |= StatusFlag::Optimal;
		case CG_ITERATIONS_EXCEED_MAXITS:
		case CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS:
			status |= StatusFlag::SubOptimal; break;
		case CG_SLOPE_ALWAYS_NEGATIVE:
			status |= StatusFlag::InsolubleModel; break;
		case CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION:
		case CG_WOLFE_CONDITIONS_NOT_SATISFIED:
		case CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES:
		case CG_NO_COST_OR_GRADIENT_IMPROVEMENT:
			status |= StatusFlag::IllConditioned; break;
		case CG_OUT_OF_MEMORY:
			status |= StatusFlag::OutOfMemory; break;
		default:
			status |= StatusFlag::Unknown; break;
		}

		return status;
	}
	void resetStatus() { status = 0; }

	void setEvaluation(CgEval value, CgEval grad) {
		model->value = value;
		model->grad = grad;
	}
	void setEvaluation(CgEval value, CgEval grad, CgEval2 valgrad) {
		setEvaluation(value, grad);
		model->valgrad = valgrad;
	}

	Index getDimension() const { return sCast<Index>(x.size()); }
	void setDimension(Index size) { x.resize(size); }

	Decision getValue(Index i) const { return x[i]; }
	void setValue(Index i, Decision value) { x[i] = value; }
	const Vec<Decision>& getValues() const { return x; }
	void setValues(const Vec<Decision>& values) { x = values; }

	void setOutput(bool shouldEnable = true) { model->Parm->PrintStatus = shouldEnable; }
	void setStatistics(bool shouldEnable = false) { model->Parm->PrintStat = shouldEnable; }
	#pragma endregion Method

protected:
	Model model;
	StatusFlags status = 0;
	Vec<Decision> x;
};

}

//#endif // _PLUGIN_CG_DESCENT


#endif // CN_HUST_GOAL_THIRD_PARTY_CG_DESCENT_H
