#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<bool> Vbool;
const ld INF = 1e17;
const ld TOL = 1e-10;
const ld PI = acos(-1);
const int LEN = 2e3;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
#define LO x
#define HI y

int N, R;
ld C[LEN]; int vp;
struct Info {
	int i;
	ld c;
	Info(int i_ = 0, ld c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return zero(c - x.c) ? i < x.i : c > x.c; }
};
std::vector<Info> G[LEN];
std::priority_queue<Info> Q;
ld dijkstra(const int& v, const int& g) {
	for (int i = 0; i < LEN; i++) C[i] = INF;
	Q.push(Info(v, 0));
	C[v] = 0;
	while (Q.size()) {
		Info p = Q.top(); Q.pop();
		if (p.c > C[p.i]) continue;
		for (Info& w : G[p.i]) {
			ld cost = p.c + w.c;
			if (C[w.i] > cost) {
				C[w.i] = cost;
				Q.push(Info(w.i, cost));
			}
		}
	}
	return C[g];
}