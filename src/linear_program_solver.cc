#include "linear_program_solver.h"

#include <cstdio>
#include <set>
#include <cassert>
#include <exception>

using std::ignore;
using std::invalid_argument;
using std::set;
using std::tie;
using std::tuple;
using std::vector;

bool LinearProgramSolver::Adjust()
{
    bool need_adjust = false;
    for (u32 i = 0; i < m; i++)
        if (coe[i][n + m] < 0)
            need_adjust = true;
    if (!need_adjust)
        return true;

    LinearProgramSolver aux_solver;

    u32 min_b_id = m;
    {
        vector<f64> c(n + 1);
        fill(c.begin() + 1, c.end(), 0);
        c[0] = -1;
        vector<vector<f64>> a(m);
        vector<f64> b(m);
        f64 min_b;
        for (u32 i = 0; i < m; i++)
        {
            a[i].resize(n + 1);
            a[i][0] = -1;
            for (u32 j = 0; j < n; j++)
                a[i][j + 1] = -coe[i][j];
            b[i] = coe[i][n + m];
            if (min_b_id >= m || (min_b_id < m && b[i] < min_b))
            {
                min_b_id = i;
                min_b = b[i];
            }
        }
        aux_solver.Reset(a, b, c);
    }
    aux_solver.Pivot(min_b_id, 0);

    f64 ans;
    tie(ignore, ans) = aux_solver.Solve();
    if (ans < 0)
        return false;

    set<u32> need_swap_in;

    for (u32 i = 1; i < m + n + 1; i++)
    {
        if (aux_solver.is_b[i] && !is_b[i - 1])
            need_swap_in.insert(i - 1);
    }

    for (u32 col_id : need_swap_in)
    {
        u32 row_id = m;
        {
            for (u32 i = 0; i < m; i++)
                if (coe[i][col_id] != 0 && !aux_solver.is_b[id[i] + 1])
                {
                    row_id = i;
                    break;
                }
        }
        assert(row_id < m);
        Pivot(row_id, col_id);
    }

    return true;
}

void LinearProgramSolver::Log() const
{
    printf("z = ");
    for (u32 i = 0; i < n + m; i++)
        if (!is_b[i])
            printf("%lf x[%d] + ", c[i], i);
    printf("%lf\n", c[n + m]);

    for (u32 i = 0; i < m; i++)
    {
        printf("x[%d] = ", id[i]);
        for (u32 j = 0; j < n + m; j++)
            if (!is_b[j])
                printf("%lf x[%d] + ", coe[i][j], j);
        printf("%lf\n", coe[i][n + m]);
    }

    printf("B = { ");
    for (u32 i = 0; i < n + m; i++)
        if (is_b[i])
            printf("%d, ", i);
    printf("}\n");
}

void LinearProgramSolver::Pivot(u32 row_id, u32 col_id)
{
    u32 old_id_row_id = id[row_id];
    f64 divider = -coe[row_id][col_id];

    id[row_id] = col_id;
    coe[row_id][col_id] = 0;
    coe[row_id][old_id_row_id] = -1;
    for (u32 i = 0; i <= n + m; i++)
        coe[row_id][i] /= divider;

    for (u32 i = 0; i < m; i++)
    {
        if (i == row_id)
            continue;
        for (u32 j = 0; j <= n + m; j++)
            coe[i][j] += coe[i][col_id] * coe[row_id][j];
        coe[i][col_id] = 0;
    }

    for (u32 i = 0; i <= n + m; i++)
        c[i] += c[col_id] * coe[row_id][i];
    c[col_id] = 0;

    is_b[old_id_row_id] = false;
    is_b[col_id] = true;
}

void LinearProgramSolver::Reset(const vector<vector<f64>> &a,
                                const vector<f64> &b, const vector<f64> &c)
{
    m = a.size();
    if (m == 0)
        throw invalid_argument("a");
    n = a.front().size();
    if (n == 0)
        throw invalid_argument("a");
    for (u32 i = 1; i < m; i++)
        if (a[i].size() != n)
            throw invalid_argument("a");
    if (b.size() != m)
        throw invalid_argument("b");
    if (c.size() != n)
        throw invalid_argument("c");

    coe.resize(m);
    id.resize(m);
    for (u32 i = 0; i < m; i++)
    {
        id[i] = i + n;
        coe[i].resize(n + m + 1);
        fill(coe[i].begin(), coe[i].end(), 0);
        for (u32 j = 0; j < n; j++)
            coe[i][j] = -a[i][j];
        coe[i][n + m] = b[i];
    }

    this->c.resize(n + m + 1);
    fill(this->c.begin(), this->c.end(), 0);
    copy(c.begin(), c.end(), this->c.begin());

    is_b.resize(n + m);
    fill(is_b.begin(), is_b.begin() + n, false);
    fill(is_b.begin() + n, is_b.end(), true);
}

tuple<LinearProgramSolver::ResultType, f64> LinearProgramSolver::Solve()
{
    using ReturnType = tuple<ResultType, f64>;

    if (!Adjust())
        return ReturnType(ResultType::kInfeasible, 0);

    while (true)
    {
        u32 col_id = n + m;
        {
            for (u32 i = 0; i < n + m; i++)
                if (!is_b[i] && c[i] > 0)
                {
                    col_id = i;
                    break;
                }
            if (col_id >= n + m)
                return ReturnType(ResultType::kSolved, c[n + m]);
        }
        u32 row_id = m;
        {
            f64 min_delta;
            for (u32 i = 0; i < m; i++)
            {
                if (coe[i][col_id] >= 0)
                    continue;
                f64 tmp = coe[i][n + m] / (-coe[i][col_id]);
                if (row_id >= m || (row_id < m && min_delta > tmp))
                {
                    min_delta = tmp;
                    row_id = i;
                }
            }
            // coe[:, col_id] >= 0
            if (row_id >= m)
                return ReturnType(ResultType::kUnbounded, 0);
        }
        Pivot(row_id, col_id);
    }
}
