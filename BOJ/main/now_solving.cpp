#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
#include <deque>
#include <random>
#include <array>
#include <tuple>
#include <complex>
#include <set>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-7;
const ld PI = acos(-1);
const int LEN = 1 << 5;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
#define LO x
#define HI y
#define LINE 1
#define CIRCLE 2
#define STRONG 0
#define WEAK 1
#define ABS 0
#define REL 1

int N, M, K, T, Q;
struct Pos {
	ld x, y;
	int i, j;
	Pos(ld x_ = 0, ld y_ = 0, int i_ = -1, int j_ = -1) : x(x_), y(y_), i(i_), j(j_) {  }
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return *this < p || *this == p; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y, i }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y, i }; }
	Pos operator * (const ld& n) const { return { x * n, y * n, i }; }
	Pos operator / (const ld& n) const { return { x / n, y / n, i }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y, i }; }
	Pos operator ~ () const { return { -y, x, i }; }
	Pos operator ! () const { return { y, x, i }; }
	ld xy() const { return x * y; }
	Pos rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2l(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} V[LEN * LEN * 4 * 4]; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
typedef std::set<Pos> MapPos;
MapPos S;
std::istream& operator >> (std::istream& is, Polygon& P) { for (Pos& p : P) is >> p.x >> p.y; return is; }
Polygon H[LEN];
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool cmpr(const Pos& p, const Pos& q) {
	bool f1 = O < p;
	bool f2 = O < q;
	if (f1 != f2) return f1;
	ld tq = p / q;
	//return zero(tq) ? p.Euc() < q.Euc() : tq > 0;
	return tq > 0;
}
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
ld dist(const Pos& d1, const Pos& d2, const Pos& q, const bool& f = REL) {
	if (f == REL) return cross(d1, d2, q) / (d1 - d2).mag();
	if (sign(dot(d1, d2, q)) <= 0 && sign(dot(d2, d1, q)) <= 0)
		return std::abs(cross(d1, d2, q)) / (d1 - d2).mag();
	return std::min((d1 - q).mag(), (d2 - q).mag());
}
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
ld area(const Polygon& H) {
	ld a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
bool norm(Polygon& H) {
	ld a = area(H);
	if (a < 0) { std::reverse(H.begin(), H.end()); return 1; }
	return 0;
}
ld rad(const Pos& d1, const Pos& d2, const Pos& d3) { return rad(d1 - d2, d3 - d2); }
int inner_check(const Polygon& H, const Pos& p) {//concave
	int cnt = 0, sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (zero(cur.y - nxt.y)) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y - TOL < p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return (cnt & 1) * 2;
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
Polygon polygon_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
	Polygon qs;
	int n = ps.size();
	for (int i = 0; i < n; i++) {
		Pos p1 = ps[i], p2 = ps[(i + 1) % n];
		int d1 = ccw(b1, b2, p1), d2 = ccw(b1, b2, p2);
		if (d1 >= 0) qs.push_back(p1);
		if (d1 * d2 < 0) qs.push_back(intersection(p1, p2, b1, b2));
	}
	return qs;
}
Polygon sutherland_hodgman(const Polygon& C, const Polygon& clip) {
	int sz = clip.size();
	std::vector<Pos> ret = C;
	for (int i = 0; i < sz; i++) {
		Pos b1 = clip[i], b2 = clip[(i + 1) % sz];
		ret = polygon_cut(ret, b1, b2);
	}
	//if (zero(area(ret))) {
	//	std::cout << "DEBUG::\n";
	//	for (Pos& p : ret) std::cout << p << "\n";
	//	std::cout << "DEBUG::\n";
	//}
	return ret;
}
Polygon circle_seg_intersection(const Pos& o, const ld& r, const Pos& p1, const Pos& p2) {
	ld d = dist(p1, p2, o);
	if (std::abs(d) > r) return {};
	Pos vec = p2 - p1;
	Pos m = intersection(p1, p2, o, o + ~vec);
	ld distance = vec.mag();
	ld ratio = sqrt(r * r - d * d);
	Pos m1 = m - vec * ratio / distance;
	Pos m2 = m + vec * ratio / distance;
	if (dot(p1, p2, m1, m2) < 0) std::swap(m1, m2);
	Polygon ret;
	if (on_seg_strong(p1, p2, m1)) ret.push_back(m1);
	if (on_seg_strong(p1, p2, m2)) ret.push_back(m2);
	return ret;//p1->p2
}
ld C[LEN * LEN * 4 * 4]; int vp;
struct Info {
	int i;
	ld c;
	Info(int i_ = 0, ld c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return c > x.c; }
};
std::vector<Info> G[LEN * LEN * 4 * 4];
ld dijkstra(const int& v, const int& g) {
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
	return C[g];
}
Polygon RV[LEN << 2];//revolve
ld get_theta(const Pos& d1, const Pos& d2, const ld& r) { return asin(r / (d1 - d2).mag()); }
Pos mid(const Pos& d1, const Pos& d2) { return (d1 + d2) * .5; }
bool between(const Pos& d1, const Pos& d2, const Pos& q) { return sign(dot(d1, d2, q)) <= 0 && sign(dot(d2, d1, q)) <= 0; }
bool close(const Pos& d, const Pos& q, const ld& r) { return (d - q).mag() < r; }
bool close(const Pos& d1, const Pos& d2, const Pos& q, const ld& r) { return dist(d1, d2, q, ABS) < r; }
Pos rotate(const Pos& p, const Pos& pv, const ld& t, const int& i) {
	Pos v = p - pv;
	ld rat = cos(t);
	Pos q = v.rot(t) * rat; q += pv; q.i = i;
	return q;
}
bool circle_is_ok(const Pos& c, const ld& r) {
	for (int j = 0; j < N; j++) {
		const Polygon& P = H[j];
		for (int i = 0; i < 4; i++) {
			if (sign(dist(P[i], P[(i + 1) % 4], c, ABS) - r) < 0) return 0;
		}
		if (inner_check(P, c)) return 0;
	}
	return 1;
}
bool connectable(const Pos& s, const Pos& e, const ld& r) {
	if (s == e) return 1;
	Pos v = ~(e - s).unit() * r;
	Polygon B = { s + v, s - v, e - v, e + v };
	for (int i = 0; i < N; i++) {
		ld a = area(sutherland_hodgman(H[i], B));
		if (!zero(a)) return 0;
	}
	return 1;
}
void connect_node(const int& n1, const int& n2, const ld& r) {
	Pos d1 = V[n1], d2 = V[n2];
	if (connectable(d1, d2, r)) {
		G[n1].push_back({ n2, (d1 - d2).mag() });
		G[n2].push_back({ n1, (d1 - d2).mag() });
	}
	return;
}
void connect_seg(const ld& r) {
	for (int i = 0; i < vp; i++)
		for (int j = i + 1; j < vp; j++)
			connect_node(i, j, r);
	return;
}
void connect_arc(const ld& r) {
	for (int n = 0; n < N; n++) {
		const Polygon& P = H[n];
		int psz = P.size();
		for (int i = 0; i < psz; i++) {
			int x = n * 4 + i;
			int sz = RV[x].size();
			for (int j = 0; j < sz; j++) RV[x][j] -= P[i];
			std::sort(RV[x].begin(), RV[x].end(), cmpr);
			for (int j = 0; j < sz; j++) RV[x][j] += P[i];
			for (int j = 0; j < sz; j++) {
				Pos lo = RV[x][j], hi = RV[x][(j + 1) % sz];
				if (lo == hi) {
					G[lo.j].push_back({ hi.j, 0 });
					G[hi.j].push_back({ lo.j, 0 });
					continue;
				}
				bool f0 = 1;
				for (int q = 0; q < N; q++) {
					const Polygon& Q = H[q];
					int qsz = Q.size();
					for (int k = 0; k < qsz; k++) {
						if (n != q || i != k) {
							bool f1 = inside(hi, P[i], lo, Q[k]);
							bool f2 = sign((P[i] - Q[k]).mag() - r * 2) < 0;
							bool f3 = sign((lo - Q[k]).mag() - r) < 0
								|| sign((hi - Q[k]).mag() - r) < 0;
							if ((f1 && f2) || f3) {
								f0 = 0;
								break;
							}
						}
						const Pos& q0 = Q[(k - 1 + psz) % psz], & q1 = Q[k];
						if (q0 == P[i] || q1 == P[i]) continue;
						Polygon inx = circle_seg_intersection(P[i], r * 2, q0, q1);
						for (const Pos& p : inx) {
							//if (inside(hi, P[i], lo, p, WEAK)) {
							if (inside(hi, P[i], lo, p, STRONG)) {
								ld d = dist(q0, q1, P[i]);
								if (d < r * 2) {
									f0 = 0;
									break;
								}
							}
						}
						if (!f0) break;
					}
					if (!f0) break;
				}
				if (f0) {
					ld t = norm(rad(lo, P[i], hi));
					ld rd = r * t;
					G[lo.j].push_back({ hi.j, rd });
					G[hi.j].push_back({ lo.j, rd });
				}
			}
		}
	}
	return;
}
void pos_init(const Pos& s, const Pos& e, const ld& r) {
	S.clear();
	vp = 0;
	V[vp++] = s;
	V[vp++] = e;
	Pos p0, p1, p2, p3;
	for (int i = 0; i < N; i++) {//tangent from S || E
		const Polygon& P = H[i];
		for (int j = 0; j < 4; j++) {
			ld t1 = get_theta(s, P[j], r);
			ld t2 = get_theta(e, P[j], r);
			p0 = rotate(P[j], s, t1, j);
			p1 = rotate(P[j], s, -t1, j);
			p2 = rotate(P[j], e, t2, j);
			p3 = rotate(P[j], e, -t2, j);
			for (const Pos& p : { p0, p1, p2, p3 }) {
				if (!S.count(p) && circle_is_ok(p, r)) {
					V[vp++] = p; S.insert(p);
				}
			}
		}
	}
	for (int i = 0; i < N; i++) {
		const Polygon& P = H[i];
		for (int k = 0; k < 4; k++) {
			int k1 = (k + 1) % 4;
			Pos v = ~(P[k] - P[k1]).unit();
			p0 = P[k] + v * r;
			p1 = P[k1] + v * r;
			p2 = P[k] - v * r;
			p3 = P[k1] - v * r;
			for (const Pos& p : { p0, p1, p2, p3 }) {
				if (!S.count(p) && circle_is_ok(p, r)) {
					V[vp++] = p; S.insert(p);
				}
			}
			for (int j = i + 1; j < N; j++) {
				const Polygon& Q = H[j];
				for (int l = 0; l < 4; l++) {
					Pos v = ~(P[k] - Q[l]).unit();
					p0 = P[k] + v * r;
					p1 = Q[l] + v * r;
					p2 = P[k] - v * r;
					p3 = Q[l] - v * r;
					for (const Pos& p : { p0, p1, p2, p3 }) {
						if (!S.count(p) && circle_is_ok(p, r)) {
							V[vp++] = p; S.insert(p);
						}
					}
					ld d = (P[k] - Q[l]).mag();
					if (d > r * 2) {//tangent from m
						Pos m = mid(P[k], Q[l]);
						ld t = get_theta(m, P[k], r);
						p0 = rotate(P[k], m, t, k);
						p1 = rotate(Q[l], m, t + PI, l);
						p2 = rotate(P[k], m, -t, k);
						p3 = rotate(Q[l], m, -t - PI, l);
						for (const Pos& p : { p0, p1, p2, p3 }) {
							if (!S.count(p) && circle_is_ok(p, r)) {
								V[vp++] = p; S.insert(p);
							}
						}
					}
				}
			}
		}
	}
	//std::cout << "vp:: " << vp << "\n";
	for (int i = 0; i < vp; i++) G[i].clear();
	for (int i = 0; i < LEN * 4; i++) RV[i].clear();
	for (int j = 2; j < vp; j++) {
		V[j].j = j;
		for (int i = 0; i < N; i++) {
			const Polygon& P = H[i];
			for (int k = 0; k < 4; k++) {
				ld d = (P[k] - V[j]).mag();
				if (zero(d - r)) RV[i * 4 + k].push_back(V[j]);
			}
		}
	}
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	int r;
	std::cin >> r >> N;
	Pos s, e; std::cin >> s >> e; s.i = -1; e.i = -2;
	for (int i = 0; i < N; i++) {
		Pos dl, ur; std::cin >> dl >> ur;
		H[i] = { dl, Pos(ur.x, dl.y), ur, Pos(dl.x, ur.y) };
	}
	pos_init(s, e, r);
	connect_seg(r);
	connect_arc(r);
	ld d = dijkstra(0, 1);
	if (d > 1e16) std::cout << "no solution\n";
	else std::cout << d << "\n";
	return;
}
int main() { solve(); return 0; }//boj3607
//boj30123 27712 3607
