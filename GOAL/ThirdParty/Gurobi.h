////////////////////////////////
/// usage : 1.	basic utilities for gurobi solver.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_THIRD_PARTY_GUROBI_H
#define CN_HUST_GOAL_THIRD_PARTY_GUROBI_H


#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"
#include "GOAL/System/System.h"

#if _PLUGIN_GUROBI

#include <functional>
#include <cstring>
#include <cmath>

#include "gurobi_c++.h"


namespace goal {

class MpSolverGurobi {
public:
	#pragma region Constant
	enum VariableType {
		Bool = GRB_BINARY,
		Integer = GRB_INTEGER,
		Real = GRB_CONTINUOUS,
		SemiInt = GRB_SEMIINT, // variables can take any integer value between the specified lower and upper bounds, or a value of zero.
		SemiReal = GRB_SEMICONT // variables can take any value between the specified lower and upper bounds, or a value of zero. 
	};

	enum OptimaDirection {
		Minimize = GRB_MINIMIZE,
		Maximize = GRB_MAXIMIZE,
		Default = Minimize
	};

	// status for the most recent optimization.
	enum StatusFlag {
		Optimal = 1,                     // GRB_OPTIMAL
		SubOptimal = 1 << 1,             // GRB_SUBOPTIMAL || (SolutionCount > 0)
		FeasibleMask = Optimal | SubOptimal,

		InsolubleModel = 1 << 2,         // GRB_INFEASIBLE || GRB_INF_OR_UNBD || GRB_UNBOUNDED
		InsolubleCutoff = 1 << 3,        // GRB_CUTOFF
		ModelErrorMask = InsolubleModel | InsolubleCutoff,

		TimeLimitExceeded = 1 << 4,      // GRB_TIME_LIMIT
		SolutionLimitExceeded = 1 << 5,  // GRB_SOLUTION_LIMIT
		NodeLimitExceeded = 1 << 6,      // GRB_NODE_LIMIT
		IterationLimitExceeded = 1 << 7, // GRB_ITERATION_LIMIT
		LimitReachedMask = TimeLimitExceeded | SolutionLimitExceeded | NodeLimitExceeded | IterationLimitExceeded,

		OutOfMemory = 1 << 8,            // OUT_OF_MEMORY
		Unknown = 1 << 9,                // any other status code or error code.
		SolverErrorMask = OutOfMemory | Unknown,
	};
	using StatusFlags = int;

	static constexpr StatusFlags SolverSpecificStatusBeginBit = 16;
	static constexpr StatusFlags GurobiErrorCodeOffset = 10000;

	enum PresolveLevel {
		Auto = GRB_PRESOLVE_AUTO,
		Off = GRB_PRESOLVE_OFF,
		Conservative = GRB_PRESOLVE_CONSERVATIVE,
		Aggressive = GRB_PRESOLVE_AGGRESSIVE,
	};

	enum MipFocusMode {
		BalancedFocus = GRB_MIPFOCUS_BALANCED,              // balance between modes.
		ImproveFeasibleSolution = GRB_MIPFOCUS_FEASIBILITY, // find sub-optima as fast as possible.
		ProveOptimality = GRB_MIPFOCUS_OPTIMALITY,          // reach and prove the optima is mandatory.
		ImproveBound = GRB_MIPFOCUS_BESTBOUND,              // improve the bound as much as possible.
	};

	enum IisMethod {
		FastIis = 0,
		SmallIis = 1,
		IgnoreBoundIis = 2, // ignores the bound constraints.
		LpRelaxIis = 3,     // the IIS for the LP relaxation of a MIP model if the relaxation is infeasible, even though the result may not be minimal when integrality constraints are included. 
		DefaultIisMethod = -1
	};

	enum QuadConvexity {
		Convex = 0,
		ConvexAfterPresolve = 1,
		NonConvex = 2,
		DefaultQuadConvexity = -1
	};

	static constexpr int MaxInt = GRB_MAXINT;
	static constexpr double MaxReal = GRB_INFINITY;

	struct MaxThreadNum { static constexpr int Auto = 0; };
	#pragma endregion Constant

	#pragma region Type
	using Model = GRBModel;
	using Decision = GRBVar;
	using Constraint = GRBConstr;
	using QuadConstraint = GRBQConstr;
	using GeneralConstraint = GRBGenConstr;
	using LinearExpr = GRBLinExpr;
	using QuadExpr = GRBQuadExpr;
	using Inequality = GRBTempConstr;

	class MpEvent;
	// on finding an MIP solution during optimization.
	using OnMipSln = std::function<void(MpEvent&)>;
	using OnMip = OnMipSln;

	class MpEvent : public GRBCallback {
	protected:
		friend MpSolverGurobi;

		OnMip onMip;
		OnMipSln onMipSln;

		MpEvent() {}
		MpEvent(OnMipSln onMipSolutionFound) : onMipSln(onMipSolutionFound) {}

	public:
		using GRBCallback::addCut;
		using GRBCallback::addLazy;
		void stop() { abort(); }

		double getValue(const Decision& var) { return getSolution(var); }
		double getValue(const LinearExpr& expr) {
			double value = expr.getConstant();
			int itemNum = sCast<int>(expr.size());
			for (int i = 0; i < itemNum; ++i) {
				value += (expr.getCoeff(i) * getValue(expr.getVar(i)));
			}
			return value; // OPTIMIZE[szx][9]: try `return expr.getValue();`?
		}
		long long getIntValue(const Decision& var) { return std::llround(getValue(var)); }
		long long getIntValue(const LinearExpr& expr) { return std::llround(getValue(expr)); }
		bool isTrue(const Decision& var) { return (getValue(var) > 0.5); }
		double getRelaxedValue(const Decision& var) { return getNodeRel(var); }

		void setValue(Decision& var, double value) { setSolution(var, value); }
		//using GRBCallback::useSolution;

		double getObj() { return getDoubleInfo(GRB_CB_MIPSOL_OBJ); }
		double getBestObj() {
			switch (where) {
			case GRB_CB_MIP:
				return getDoubleInfo(GRB_CB_MIP_OBJBST);
			case GRB_CB_MIPNODE:
				return getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
			case GRB_CB_MIPSOL:
				return getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
			default:
				return 0;
			}
		}
		double getBestBound() {
			switch (where) {
			case GRB_CB_MIP:
				return getDoubleInfo(GRB_CB_MIP_OBJBND);
			case GRB_CB_MIPNODE:
				return getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
			case GRB_CB_MIPSOL:
				return getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
			default:
				return 0;
			}
		}

		void callback() {
			if (where == GRB_CB_MIP) {
				if (onMip) { onMip(*this); }
			} else if (where == GRB_CB_MIPSOL) {
				if (onMipSln) { onMipSln(*this); }
			}
		}
	};
	#pragma endregion Type

	#pragma region Constructor
	MpSolverGurobi() : model(getEnv()) {}

	void loadModel(const Str& inputPath) { model.read(inputPath); }
	void saveModel(const Str& outputPath = "model.lp") {
		model.update();
		model.write(outputPath);
	}
	#pragma endregion Constructor

	#pragma region Method
	StatusFlags optimize() {
		resetStatus();
		try {
			//model.update();
			model.optimize();
		} catch (GRBException& e) {
			if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
				status |= StatusFlag::OutOfMemory;
			} else {
				status |= StatusFlag::Unknown;
			}
			setSolverSpecificStatus(e.getErrorCode() - GurobiErrorCodeOffset);
		}
		return getStatus();
	}

	void computeIIS(const Str& outputPath = "iis.ilp") {
		try {
			model.computeIIS();
			model.write(outputPath);
		} catch (GRBException& e) {
			status |= StatusFlag::Unknown;
			setSolverSpecificStatus(e.getErrorCode() - GurobiErrorCodeOffset);
		}
	}

	Decision addVar(VariableType type, double lb = 0, double ub = 1, double objCoef = 0, const Str& name = "") {
		return model.addVar(lb, ub, objCoef, sCast<char>(type), name);
	}

	Constraint addConstraint(const Inequality& ieq, const Str& name = "") { return model.addConstr(ieq, name); }
	QuadConstraint addQuadConstraint(const Inequality& ieq, const Str& name = "") { return model.addQConstr(ieq, name); }
	GeneralConstraint addAbsConstraint(Decision absVar, Decision freeVar, const Str& name = "") {
		return model.addGenConstrAbs(absVar, freeVar, name);
	}
	GeneralConstraint addIndicatorConstraint(Decision boolVar, bool value, const Inequality& ieq, const Str& name = "") {
		return model.addGenConstrIndicator(boolVar, value, ieq, name);
	}
	void removeConstraint(Constraint constraint) { model.remove(constraint); }
	int getConstraintCount() const { return model.get(GRB_IntAttr_NumConstrs); }

	void setObjective(const LinearExpr& expr, OptimaDirection optimaDirection) {
		model.setObjective(expr, optimaDirection);
	}
	void setOptimaDirection(OptimaDirection optimaDirection = OptimaDirection::Default) {
		model.set(GRB_IntAttr_ModelSense, optimaDirection);
	}

	static GRBEnv& getEnv() {
		thread_local static bool init = true;
		if (init) {
			os::StdIo::disableStdOut();
			globalEnv.start();
			//std::memset(&NoVar, 0, sizeof(Decision)); // static variable is zero-initialized by default.
			init = false;
			os::StdIo::restoreStdOut();
		}
		return globalEnv;
	}

	StatusFlags getStatus() {
		status &= StatusFlag::SolverErrorMask;

		if (model.get(GRB_IntAttr_SolCount) > 0) { status |= StatusFlag::SubOptimal; }

		switch (model.get(GRB_IntAttr_Status)) {
		case GRB_OPTIMAL:
			status |= StatusFlag::Optimal;
		case GRB_SUBOPTIMAL:
			status |= StatusFlag::SubOptimal; break;
		case GRB_LOADED:
		case GRB_INPROGRESS:
			break;
		case GRB_ITERATION_LIMIT:
			status |= StatusFlag::IterationLimitExceeded; break;
		case GRB_NODE_LIMIT:
			status |= StatusFlag::NodeLimitExceeded; break;
		case GRB_TIME_LIMIT:
			status |= StatusFlag::TimeLimitExceeded; break;
		case GRB_SOLUTION_LIMIT:
			status |= StatusFlag::SolutionLimitExceeded; break;
		case GRB_INFEASIBLE:
		case GRB_INF_OR_UNBD:
		case GRB_UNBOUNDED:
			status |= StatusFlag::InsolubleModel; break;
		case GRB_CUTOFF:
			status |= StatusFlag::InsolubleCutoff; break;
		default:
			status |= StatusFlag::Unknown; break;
		}

		return status;
	}
	StatusFlags getSolverSpecificStatus() { return status >> SolverSpecificStatusBeginBit; }
	void setSolverSpecificStatus(StatusFlags solverSpecificStatus) { status |= (solverSpecificStatus << SolverSpecificStatusBeginBit); }
	void resetStatus() { status = 0; }

	double getObjValue() const { return model.get(GRB_DoubleAttr_ObjVal); }
	long long getIntObjValue() const { return std::llround(getObjValue()); }
	double getObjBound() const { return model.get(GRB_DoubleAttr_ObjBound); }
	double getObjBoundReal() const { return model.get(GRB_DoubleAttr_ObjBoundC); }

	static double getValue(const LinearExpr& expr) { return expr.getValue(); }
	static double getValue(const Decision& var) { return var.get(GRB_DoubleAttr_X); }
	static long long getIntValue(const LinearExpr& expr) { return std::llround(getValue(expr)); }
	static long long getIntValue(const Decision& var) { return std::llround(getValue(var)); }

	static bool isTrue(double value) { return (value > 0.5); }
	static bool isTrue(LinearExpr expr) { return isTrue(getValue(expr)); }
	static bool isTrue(Decision var) { return isTrue(getValue(var)); }

	// make sure `Decision` variables have been set to 0 before using this method.
	static bool initialized(const Decision& var) { return std::memcmp(&var, &NoVar, sizeof(Decision)) != 0; }

	// the methods in MpSolver is invalid within the callback, only use the ones in MpEvent instead.
	void setMipSlnEvent(OnMipSln onMipSln, bool enbaleLazyConstraints = true) {
		if (enbaleLazyConstraints) { model.set(GRB_IntParam_LazyConstraints, 1); }
		mpEvent.onMipSln = onMipSln;
		model.setCallback(&mpEvent);
	}
	void setMipEvent(OnMipSln onMip) {
		mpEvent.onMip = onMip;
		model.setCallback(&mpEvent);
	}

	void setTimeLimit(double secTimeout) {
		if (secTimeout < 0) { secTimeout = 0; }
		model.set(GRB_DoubleParam_TimeLimit, secTimeout);
	}
	void setBestBoundToStop(double bound) { model.set(GRB_DoubleParam_BestBdStop, bound); }
	void setBestObjToStop(double bound) { model.set(GRB_DoubleParam_BestObjStop, bound); }
	void setOutput(bool shouldEnable = true) { model.set(GRB_IntParam_OutputFlag, shouldEnable); }
	void setMaxThread(int threadNum = MaxThreadNum::Auto) { model.set(GRB_IntParam_Threads, threadNum); }
	void setSeed(int seed) { model.set(GRB_IntParam_Seed, (seed & (std::numeric_limits<int>::max)())); }

	void setPresolveLevel(PresolveLevel presolveLevel) { model.set(GRB_IntParam_Presolve, presolveLevel); }
	void setMipFocus(MipFocusMode mode) { model.set(GRB_IntParam_MIPFocus, mode); }

	double getLB(Decision& var) const { return var.get(GRB_DoubleAttr_LB); }
	double getUB(Decision& var) const { return var.get(GRB_DoubleAttr_UB); }
	void setLB(Decision& var, double value) { var.set(GRB_DoubleAttr_LB, value); }
	void setUB(Decision& var, double value) { var.set(GRB_DoubleAttr_UB, value); }
	void setBound(Decision& var, double lb, double ub) { setLB(var, lb); setUB(var, ub); }
	void setType(Decision& var, VariableType vt) { var.set(GRB_CharAttr_VType, vt); }

	// [Tune] use the given value as the initial solution in MIP.
	void setInitValue(Decision& var, double value) { var.set(GRB_DoubleAttr_Start, value); }
	// [Tune] guide the solver to prefer certain value on certain variable.
	void setHintValue(Decision& var, double value) { var.set(GRB_DoubleAttr_VarHintVal, value); }
	void setHintValue(Decision& var, double value, int confidence) {
		setHintValue(var, value);
		var.set(GRB_IntAttr_VarHintPri, confidence);
	}
	// [Tune] guide the solver to prefer certain variable for branching.
	void setBranchPriority(Decision& var, int priority) { var.set(GRB_IntAttr_BranchPriority, priority); }

	void setIisMethod(IisMethod iisMethod) { model.set(GRB_IntParam_IISMethod, iisMethod); }

	void setOptimalityTolerance(double value) { model.set(GRB_DoubleParam_OptimalityTol, value); }
	void setFeasibilityTolerance(double value) { model.set(GRB_DoubleParam_FeasibilityTol, value); }
	// tolerance for integer variable type constraints.
	void setIntegerTolerance(double value) { model.set(GRB_DoubleParam_IntFeasTol, value); }

	void setQuadConvexity(QuadConvexity quadConvexity) { model.set(GRB_IntParam_NonConvex, quadConvexity); }
	#pragma endregion Method

	static const Decision NoVar;

protected:
	thread_local static GRBEnv globalEnv;

	Model model;
	MpEvent mpEvent;
	StatusFlags status = 0;
};

}

#endif // _PLUGIN_GUROBI


#endif // CN_HUST_GOAL_THIRD_PARTY_GUROBI_H
