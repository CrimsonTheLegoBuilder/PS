#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
typedef long long ll;
typedef double ld;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-9;
const int V_LEN = 1 << 10;
const int G_LEN = 1 << 19;
const int F_LEN = 1 << 11;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ll sq(const ll& x) { return x * x; }

ld C[G_LEN];
int vp;
struct Info {
	int i;
	ld c;
	Info(int i_ = 0, ld c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return c > x.c; }
};
std::vector<Info> G[G_LEN];
void dijkstra(const int& v) {
	std::priority_queue<Info> PQ;
	for (int i = 0; i < vp; i++) C[i] = INF;
	PQ.push(Info(v, 0));
	C[v] = 0;
	while (PQ.size()) {
		Info p = PQ.top(); PQ.pop();
		if (p.c > C[p.i]) continue;
		for (Info& w : G[p.i]) {
			ld cost = p.c + w.c;
			if (C[w.i] > cost) {
				C[w.i] = cost;
				PQ.push(Info(w.i, cost));
			}
		}
	}
	return;
}
int T, N, R, M, Q;
ld L[V_LEN], lo, hi;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return sqrt(Euc()); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} S[V_LEN], E[V_LEN], qry[V_LEN], FS[F_LEN], FE[F_LEN]; ld FC[F_LEN]; int f;
const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
ld intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2, const bool& f = 0) {
	ll tq = (q2 - q1) / (p2 - p1);
	if (tq == 0) {
		if (q1 == p1 || q2 == p1) return 0;
		if (q1 == p2 || q2 == p2) return 1;
		return -1;
	}
	return ((q2 - q1) / (q1 - p1)) / (ld)tq;
}
inline bool circle_inner_check(const Pos& q, const ll& r, const Pos& p) { return sq(r) >= (p - q).Euc(); }
bool circle_line_intersections(const Pos& q, const ll& r, const Pos& s, const Pos& e, ld& lo, ld& hi) {
	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
	Pos vec = e - s;
	Pos OM = s - q;
	ll a = vec.Euc();
	ll b = vec * OM;
	ll c = OM.Euc() - sq(r);
	ll J = b * b - a * c;
	if (J < 0) return 0;
	ld det = sqrt(J);
	lo = (-b - det) / a;
	hi = (-b + det) / a;
	if (!circle_inner_check(q, r, s) && !circle_inner_check(q, r, e))
		if ((lo < 0 && hi < 0) || (lo > 1 && hi > 1)) return 0;
	return 1;
}
struct Event {
	int i;
	ld x;
	Event(int i_ = 0, ld x_ = 0) : i(i_), x(x_) {}
	bool operator < (const Event& e) const { return x < e.x; }
};
std::vector<Event> X[V_LEN];
ld ans[V_LEN];
bool V[V_LEN];
void test() {
	for (int i = 0; i < vp; i++) G[i].clear();
	for (int i = 0; i < N; i++) X[i].clear();
	memset(V, 0, sizeof V);
	vp = 1; f = 0;
	std::cin >> N >> R;
	for (int i = 0; i < N; i++) {
		std::cin >> S[i] >> E[i] >> M;
		L[i] = (S[i] - E[i]).mag();
		if (!M) continue;
		for (int j = 0; j < M; j++) {
			ld x; std::cin >> x;
			X[i].push_back(Event(vp, x));
			G[0].push_back(Info(vp, 0));
			FS[f] = S[i]; FE[f] = E[i]; FC[f] = x;
			vp++; f++;
		}
	}
	for (int i = 0; i < N; i++) {
		const Pos& si = S[i], & ei = E[i];
		for (int j = i + 1; j < N; j++) {
			const Pos& sj = S[j], & ej = E[j];
			if (!intersect(si, ei, sj, ej)) continue;
			ld x;
			x = intersection(si, ei, sj, ej); assert(x > -.5);
			X[i].push_back(Event(vp, x));
			x = intersection(sj, ej, si, ei); assert(x > -.5);
			X[j].push_back(Event(vp, x));
			vp++;
		}
	}
	for (int i = 0; i < N; i++) {
		ld l = (S[i] - E[i]).mag();
		std::sort(X[i].begin(), X[i].end());
		int sz = X[i].size();
		for (int j = 0; j < sz - 1; j++) {
			ld d = X[i][j + 1].x - X[i][j].x;
			G[X[i][j].i].push_back(Info(X[i][j + 1].i, l * d));
			G[X[i][j + 1].i].push_back(Info(X[i][j].i, l * d));
		}
	}

	dijkstra(0);

	std::cin >> Q;
	for (int k = 1; k <= Q; k++) {
		std::cin >> qry[k]; V[k] = 0; ans[k] = INF;
		for (int m = 0; m < f; m++) {
			if (!circle_line_intersections(qry[k], R, FS[m], FE[m], lo, hi)) continue;
			ld x = FC[m];
			if (lo <= x && x <= hi) { V[k] = 1; ans[k] = 0; break; }
		}
	}
	std::vector<Event> A; A.reserve(V_LEN << 1);
	for (int i = 0; i < N; i++) {
		ld & l = L[i], cs, ce, d;
		Event s = Event(0, 0), e = Event(0, 1);
		A.clear();
		for (const Event& v : X[i]) A.push_back(v);
		for (int k = 1; k <= Q; k++) {
			if (V[k]) continue;
			if (!circle_line_intersections(qry[k], R, S[i], E[i], lo, hi)) continue;
			if (-TOL < lo && lo < 1 + TOL) A.push_back(Event(-k, lo));
			if (-TOL < hi && hi < 1 + TOL) A.push_back(Event(-k, hi));
		}
		std::sort(A.begin(), A.end());
		int sz = A.size();
		cs = INF, ce = INF;
		for (int k = 0, j; k < sz; k++) {
			j = sz - k - 1;
			d = (A[k].x - s.x) * l;
			if (A[k].i < 0) ans[-A[k].i] = std::min(ans[-A[k].i], cs + d);
			else { if (cs + d > C[A[k].i]) s = A[k], cs = C[A[k].i]; }
			d = (e.x - A[j].x) * l;
			if (A[j].i < 0) ans[-A[j].i] = std::min(ans[-A[j].i], ce + d);
			else { if (ce + d > C[A[j].i]) e = A[j], ce = C[A[j].i]; }
		}
	}
	for (int k = 1; k <= Q; k++) {
		if (ans[k] > 1e16) std::cout << "-1\n";
		else std::cout << ans[k] << "\n";
	}
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(7);
	std::cin >> T; while (T--) test();
	return;
}
int main() { solve(); return 0; }//boj10732