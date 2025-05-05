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
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const int LEN = 1 << 10;
const ld TOL = 1e-7;
const ll MOD = 1e9 + 7;
const ld PI = acos(-1);
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
ld flip(ld lat) {
	if (zero(lat - PI * .5) || zero(lat + PI * .5)) return 0;
	if (zero(lat)) return PI * .5;
	if (lat > 0) return PI * .5 - lat;
	if (lat < 0) return -(PI * .5) - lat;
	return INF;
}
//ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

int N, M, T, Q;
struct Pos {
	int x, y;
	//ll x, y;
	int i, j;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1, int j_ = -1) : x(x_), y(y_), i(i_), j(j_) {}
	//Pos(ll x_ = 0, ll y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	int Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
Polygon P[10];
Pos V[LEN];
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
int collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2, const int& f = STRONG) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	if (f == WEAK) return f1 && f2;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
ll area(const Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) { if (area(H) < 0) std::reverse(H.begin(), H.end()); }
int inner_check(const Polygon& H, const Pos& p) {//concave
	int cnt = 0, sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (cur.y == nxt.y) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y <= p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return (cnt & 1) * 2;
}
ll C[LEN]; int vp;
struct Info {
	int i;
	ll c;
	Info(int i_ = 0, ll c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return c > x.c; }
};
std::vector<Info> G[LEN];
ld dijkstra(const int& v, const int& g) {
	std::priority_queue<Info> PQ;
	for (int i = 0; i < LEN; i++) C[i] = INF;
	PQ.push(Info(v, 0));
	C[v] = 0;
	while (PQ.size()) {
		Info p = PQ.top(); PQ.pop();
		if (p.c > C[p.i]) continue;
		if (p.i == g) return C[g];
		for (Info& w : G[p.i]) {
			ld cost = p.c + w.c;
			if (C[w.i] > cost) {
				C[w.i] = cost;
				PQ.push(Info(w.i, cost));
			}
		}
	}
	return C[g];
}
bool connectable(const int& i, const int& j) {
	const Pos& s = V[i], & g = V[j];
	for (int n = 0; n < N; n++) {
		const Polygon& H = P[n];
		int sz = H.size();
		for (int k = 0; k < sz; k++) {
			const Pos& k0 = H[k], k1 = H[(k + 1) % sz];
			if (intersect(s, g, k0, k1, WEAK)) return 0;
			if (on_seg_weak(s, g, k0)) return 0;
			if (on_seg_weak(s, g, k1)) return 0;
		}
	}
	if (V[i].i == V[j].i) {
		assert(V[i].i > -1);
		const Polygon& H = P[V[i].i];
		int sz = H.size();
		int u = V[i].j;
		int v = V[j].j;
		if (u == (v + 1) % sz || v == (u + 1) % sz) return 1;
		const Pos& p0 = H[(u - 1 + sz) % sz], & p1 = H[u], & p2 = H[(u + 1) % sz];
		const Pos& q0 = H[(v - 1 + sz) % sz], & q1 = H[v], & q2 = H[(v + 1) % sz];
		if (inside(p0, p1, p2, q1)) return 0;
		if (inside(q0, q1, q2, p1)) return 0;
	}
	return 1;
}
ld query(const Pos& s, const Pos& g) {
	V[0] = s; V[1] = g;
	G[0].clear();
	if (connectable(0, 1)) {
		int d = std::max(V[1].y - V[0].y, 0);
		return d;
	}
	for (int i = 2; i < vp; i++)
		if (G[i].size() && G[i].back().i == 1)
			G[i].pop_back();
	for (int i = 2; i < vp; i++) {
		if (connectable(0, i)) {
			int d = V[i].y - V[0].y;
			G[0].push_back(Info(i, std::max(d, 0)));
		}
		if (connectable(1, i)) {
			int d = V[1].y - V[i].y;
			G[i].push_back(Info(1, std::max(d, 0)));
		}
	}
	ll d = dijkstra(0, 1);
	return d;
}
bool query() {
	std::cin >> N >> Q;
	if (!N && !Q) return 0;
	for (int i = 0; i < vp; i++) G[i].clear();
	vp = 2;
	for (int i = 0; i < N; i++) {
		int v; std::cin >> v;
		P[i].resize(v);
		for (Pos& p : P[i]) { std::cin >> p; p.i = i; }
		norm(P[i]);
		int j = 0;
		for (Pos& p : P[i]) { p.j = j++; V[vp++] = p; }
	}
	for (int i = 2; i < vp; i++) {
		for (int j = i + 1; j < vp; j++) {
			if (connectable(i, j)) {
				int d = V[j].y - V[i].y;
				G[i].push_back(Info(j, std::max(d, 0)));
				d *= -1;
				G[j].push_back(Info(i, std::max(d, 0)));
			}
		}
	}
	while (Q--) {
		Pos s, g; std::cin >> s >> g;
		s.i = -1; g.i = -2;
		std::cout << query(s, g) << "\n";
	}
	return 1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	while (query());
	return;
}
int main() { solve(); return 0; }//boj22367

/*

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <cstdlib>
using namespace std;

#define ALL(v) begin(v),end(v)
namespace {

	using LL = long long;
	using pll = pair<LL, LL>;
	using segment = pair<pll, pll>;
	using vpll = vector<pll>;
	using D = long double;
	using point = pair<D, LL>;

	const LL INF = 1010101010;

	D middle(D x, D y) { return (x + y) / 2; }

	struct solver {
		vector<vpll> pols;
		vector<segment> segs;

		static bool calcx(pll p1, pll p2, LL y, D& res) {
			if (p1.second == p2.second) { return false; }
			if (p1.second > p2.second) { swap(p1, p2); }
			if (y < p1.second || y > p2.second) { return false; }
			LL dx = p2.first - p1.first, dy = p2.second - p1.second;
			res = D(y - p1.second) * dx / dy + p1.first;
			return true;
		}

		bool inanypol(D x, LL y) {
			for (const auto& pol : pols) {
				int cnt = 0;
				int m = pol.size() - 1;
				for (int i = 0; i < m; ++i) {
					D r;
					pll p1 = pol[i], p2 = pol[i + 1];
					if (p1.second == p2.second && p1.second == y) {
						LL x1, x2;
						tie(x1, x2) = minmax(p1.first, p2.first);
						D r1 = x1, r2 = x2;
						if (r1 <= x && x <= r2) {
							if (r1 < x && x < r2) { return false; }
							if (right && x == r1) { return false; }
							if (!right && x == r2) { return false; }
						}
					}
					if (max(p1.second, p2.second) <= y) { continue; }
					if (calcx(p1, p2, y, r)) {
						if (r < x) { ++cnt; }
					}
				}
				if (cnt & 1) { return true; }
			}
			return false;
		}

		void solve(int npol) {
			int q;
			cin >> q;

			vector<point> pts;
			map<point, int> pttoidx;
			vector<vector<pll>> G;
			auto regpt = [&](const point& pt) {
				int id = pttoidx.emplace(pt, pttoidx.size()).first->second;
				if (id >= (int)pts.size()) {
					pts.push_back(pt);
					G.emplace_back();
				}
				return id;
				};

			map<LL, vector<int>> ytoids;

			segs.emplace_back(pll(-INF, -INF), pll(-INF, INF));
			segs.emplace_back(pll(INF, -INF), pll(INF, INF));
			for (int i = 0; i < npol; ++i) {
				int m;
				cin >> m;
				vector<pll> pl(m + 1);
				for (int j = 0; j < m; ++j) {
					cin >> pl[j].first >> pl[j].second;
					ytoids[pl[j].second].push_back(regpt(point(pl[j].first, pl[j].second)));
				}
				pl[m] = pl[0];

				for (int j = 0; j < m; ++j) {
					pll p1 = pl[j], p2 = pl[j + 1];
					if (p1.second > p2.second) { swap(p1, p2); }
					if (p1.second != p2.second) {
						segs.emplace_back(p1, p2);
					}
					else {
						int u1 = regpt(point(p1.first, p1.second));
						int u2 = regpt(point(p2.first, p2.second));
						G[u1].emplace_back(u2, 0);
						G[u2].emplace_back(u1, 0);
					}
				}

				pols.emplace_back(move(pl));
			}

			vector<int> starts(q), goals(q);
			for (int i = 0; i < q; ++i) {
				LL sx, sy, gx, gy;
				cin >> sx >> sy >> gx >> gy;
				starts[i] = regpt(point(sx, sy));
				goals[i] = regpt(point(gx, gy));

				ytoids[sy].push_back(starts[i]);
				ytoids[gy].push_back(goals[i]);
			}

			int ns = segs.size();
			vector<vector<int>> idsonsegs(ns);
			for (int i = 0; i < ns; ++i) {
				const segment& s = segs[i];
				idsonsegs[i].push_back(regpt(point(s.first.first, s.first.second)));
			}
			for (auto& pr : ytoids) {
				LL y = pr.first;
				vector<pair<D, int>> itss;

				for (int i = 0; i < ns; ++i) {
					D x;
					if (calcx(segs[i].first, segs[i].second, y, x)) {
						itss.emplace_back(x, i);
					}
				}

				for (int id : pr.second) {
					pair<D, int> prv(-INF, 0), nxt(INF, 1);
					for (const auto& p : itss) {
						if (p.first < pts[id].first) {
							prv = max(prv, p);
						}
						if (pts[id].first < p.first) {
							nxt = min(nxt, p);
						}
					}

					if (!inanypol(middle(prv.first, pts[id].first), pts[id].second)) {
						int k = regpt(point(prv.first, y));
						G[id].emplace_back(k, 0);
						G[k].emplace_back(id, 0);
						idsonsegs[prv.second].push_back(k);
					}
					if (!inanypol(middle(pts[id].first, nxt.first), pts[id].second)) {
						int k = regpt(point(nxt.first, y));
						G[id].emplace_back(k, 0);
						G[k].emplace_back(id, 0);
						idsonsegs[nxt.second].push_back(k);
					}
				}
			}

			for (int i = 0; i < ns; ++i) {
				regpt(point(segs[i].second.first, segs[i].second.second));
				idsonsegs[i].push_back(regpt(point(segs[i].second.first, segs[i].second.second)));
				const auto& v = idsonsegs[i];
				for (size_t j = 0; j + 1 < v.size(); ++j) {
					int id1 = v[j];
					int id2 = v[j + 1];
					G[id1].emplace_back(id2, pts[id2].second - pts[id1].second);
					G[id2].emplace_back(id1, 0);
				}
			}

			priority_queue<pll> pq;
			vector<LL> dist;
			for (int i = 0; i < q; ++i) {
				pq.emplace(0, starts[i]);
				dist.assign(G.size(), 1LL << 60);
				dist[starts[i]] = 0;
				while (!pq.empty()) {
					LL d = -pq.top().first;
					int u = pq.top().second;
					pq.pop();
					if (dist[u] != d) { continue; }
					for (const pll& e : G[u]) {
						LL nd = d + e.second;
						if (dist[e.first] > nd) {
							dist[e.first] = nd;
							pq.emplace(-nd, e.first);
						}
					}
				}
				cout << dist[goals[i]] << '\n';
			}
		}
	};
}

int main() try {
	int n;
	while (cin >> n && n) {
		solver().solve(n);
	}
}
catch (const exception& ex) {
	cerr << ex.what() << endl;
	abort();
}

*/