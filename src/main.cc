
#include "linear_program_solver.h"
#include "timer.h"
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

struct edge {
    int u, v;
    int w;

    bool operator<(const edge& b) const {
        const int a1[] = {u, v, w};
        const int a2[] = {b.u, b.v, b.w};
        for (auto i = 0; i < 3; i++) {
            if (a1[i] < a2[i]) {
                return true;
            } else if (a1[i] > a2[i]) {
                return false;
            }
        }
        return false;
    }
};

constexpr int REP = 100;
constexpr int MAX_D = 100;
constexpr int INF = (int) 1e6;

struct ShortestPathBenchmark {
    vector<edge> edge_lst;
    int n;  // number of vertices
    int m;  // number of edges

    int sources[REP];
    int targets[REP];
    double result_dij[REP];
    double result_lp[REP];

    void generate(int n) {
        this->n = n;
        this->m = (int) (n * sqrt(n));

        map<edge, int> edges;
        while (edges.size() < m) {
            auto u = rand() % n;
            auto v = rand() % n;
            if (u != v) {
                edge e{u, v, 0};
                edges[e] = rand() % MAX_D;
            }
        }

        edge_lst.reserve(m);

        for (auto e:edges) {
            edge_lst.push_back({e.first.u, e.first.v, e.second});
        }

        int p = 0;
        while (p < REP) {
            sources[p] = rand() % n;
            targets[p] = rand() % n;
            if (sources[p] != targets[p]) p++;
        }
    }

    void dijkstra() {

        // build
        vector<list<edge>> adj_list(n);
        for (auto e:edge_lst) {
            adj_list[e.u].push_back({e.u, e.v, e.w});
            adj_list[e.v].push_back({e.v, e.u, e.w});
        }

        timer timer;
        for (auto r = 0; r < REP; r++) {
            auto s = sources[r];
            auto t = targets[r];

            // init
            int dis[n];
            bool inset[n];
            for (auto i = 0; i < n; i++) {
                dis[i] = INF;
                inset[i] = false;
            }
            dis[s] = 0;

            // main algorithm
            for (auto k = 0; k < n; k++) {

                int min_d = INF, min_idx = -1;
                for (auto idx = 0; idx < n; idx++) {
                    auto d = dis[idx];
                    if (!inset[idx] && d < min_d) {
                        min_d = d;
                        min_idx = idx;
                    }
                }

                if (min_idx == -1) break;

                inset[min_idx] = true;
                for (auto e:adj_list[min_idx]) {
                    auto new_dis = dis[e.u] + e.w;
                    if (new_dis < dis[e.v]) {
                        dis[e.v] = new_dis;
                    }
                }

                if (min_idx == t) break;
            }

            result_dij[r] = dis[t];
        }

        auto s = timer.get_seconds();
        cout << "dijkstra, " << n << ", " << s << endl;
    }

    void lp() {
        vector<double> b;
        vector<double> c;
        vector<vector<double>> a;
        timer timer;

        for (auto r = 0; r < REP; r++) {
            LinearProgramSolver solver;
            auto s = sources[r];
            auto t = targets[r];

            // target d_t
            c.resize(n);
			fill(c.begin(), c.end(), 0);
            c[t] = 1;

            // restrictions

            // d_s = 0
			b.clear();
            b.push_back(0);
            b.push_back(0);

			a.clear();
            vector<double> vec(n);
            vec[s] = 1;
            a.push_back(vec);
            vec[s] = -1;
            a.push_back(vec);

            // d_v <= d_u + w(u,v)
            for (auto e:edge_lst) {
                vec.resize(n);
				fill(vec.begin(), vec.end(), 0);
                vec[e.v] = 1;
                vec[e.u] = -1;
                a.push_back(vec);
                b.push_back(e.w);

                fill(vec.begin(), vec.end(), 0);
				vec[e.v] = -1;
                vec[e.u] = 1;
                a.push_back(vec);
                b.push_back(e.w);
            }
            solver.Reset(a, b, c);
			LinearProgramSolver::ResultType result_type;
			tie(result_type, result_lp[r]) = solver.Solve();
		}

        auto s = timer.get_seconds();
        cout << "linear programming, " << n << ", " << s << endl;
    }

    void run() {
        auto n = 4;
        for (auto i = 0; i < 9; i++) {
            generate(n);
            dijkstra();
            lp();
            for (auto j = 0; j < REP; j++) {
                auto delta = result_dij[j] - result_lp[j];
                if (fabs(delta) > 1e-5) {
                    cout << "err = " << delta
						<< ", result_dij = " << result_dij[j] 
						<< ", result_lp = " << result_lp[j] << "\n";
                }
            }
            n *= 2;
        }
    }
};


ShortestPathBenchmark benchmark;

int main() {
    benchmark.run();
    return 0;
}