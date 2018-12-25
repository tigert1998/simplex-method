#pragma once

#include "types.h"

#include <vector>
#include <tuple>

class LinearProgramSolver
{
  private:
	u32 n, m;

	// coe: m * (n + m + 1)
	// id: m
	// c: (n + m + 1)
	// is_b: (n + m)
	std::vector<std::vector<f64>> coe;
	std::vector<u32> id;
	std::vector<f64> c;
	std::vector<bool> is_b;

	bool Adjust();

	void Log() const;

	void Pivot(u32 row_id, u32 col_id);

  public:
	enum class ResultType
	{
		kUnbounded,
		kSolved,
		kInfeasible
	};

	void Reset(const std::vector<std::vector<f64>> &a,
			   const std::vector<f64> &b, const std::vector<f64> &c);

	std::tuple<ResultType, double> Solve();
};