
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

constexpr int REP = 1000;
constexpr int MAX_D = 100;
constexpr int INF = (int) 1e6;

struct ShortestPathBenchmark {
    vector<edge> edge_lst;
    int n;  // number of vertices
    int m;  // number of edges

    int sources[REP];
    int result_dij[REP];

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

        for (int& source : sources) {
            source = rand() % n;
        }
    }

    void dijkstra() {

        // build
        vector<list<edge>> adj_list(n);
        for (auto e:edge_lst) {
            adj_list[e.u].push_back({e.u, e.v, e.w});
            adj_list[e.v].push_back({e.v, e.u, e.w});
        }

        timer t;
        for (auto r = 0; r < REP; r++) {
            auto s = sources[r];

            // build
            int dis[n];
            bool inset[n];
            for (auto i = 0; i < n; i++) {
                dis[i] = INF;
                inset[i] = false;
            }
            dis[s] = 0;

            int max_d = 0;  // as certificate

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
                max_d = max(max_d, min_d);

                if (min_idx == -1) break;

                inset[min_idx] = true;
                for (auto e:adj_list[min_idx]) {
                    auto new_dis = dis[e.u] + e.w;
                    if (new_dis < dis[e.v]) {
                        dis[e.v] = new_dis;
                    }
                }
            }
            result_dij[r] = max_d;
        }

        auto s = t.get_seconds();
        cout << "dijkstra, " << n << ", " << s << endl;
    }

    void lp() {

    }

    void run() {
        auto n = 4;
        for (auto i = 0; i < 9; i++) {
            generate(n);
            dijkstra();
            n *= 2;
        }
    }
};


ShortestPathBenchmark benchmark;

int main() {
    benchmark.run();
    return 0;
}