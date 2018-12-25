#include "linear_program_solver.h"

LinearProgramSolver solver;

int main() {
	using std::vector;
	vector<double> b, c;
	vector<vector<double>> a;

	int n, m;
	scanf("%d%d", &n, &m);
	b.resize(m);
	c.resize(n);
	a.resize(m);
	for (int i = 0; i < n; i++) 
		scanf("%lf", &c[i]);
	for (int i = 0; i < m; i++) {
		a[i].resize(n);
		for (int j = 0; j < n; j++) 
			scanf("%lf", &a[i][j]);
		scanf("%lf", &b[i]);
	}
	solver.Reset(a, b, c);
	double ans;
	LinearProgramSolver::ResultType result_type;

	std::tie(result_type, ans) = solver.Solve();

	printf("%.6lf", ans);
	return 0;
}