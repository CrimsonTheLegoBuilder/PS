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
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-7;
const ld PI = acos(-1);
const int LEN = 1e3;
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
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
//ll gcd(ll a, ll b) {
//	while (b) {
//		ll tmp = a % b;
//		a = b;
//		b = tmp;
//	}
//	return a;
//}

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2
#define SEG 3
#define POS 4

#define STRONG 0
#define WEAK 1

//freopen("../../../input_data/triathlon_tests/triath.20", "r", stdin);
//freopen("../../../input_data/triathlon_tests/triathlon_out.txt", "w", stdout);

//Euler characteristic : v - e + f == 1
//Pick`s Theorem : A = i + b/2 - 1

//2D============================================================================//

int N, M, K;
struct Pos {
	ld x, y;
	Pos(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return *this < p || *this == p; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& n) const { return { x * n, y * n }; }
	Pos operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pos rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	bool close(const Pos& p) const { return zero((*this - p).Euc()); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
ld dist(const Pos& d1, const Pos& d2, const Pos& t, bool f = 0) {
	if (!f) return cross(d1, d2, t) / (d1 - d2).mag();
	if (sign(projection(d1, d2, d2, t)) <= 0 &&
		sign(projection(d2, d1, d1, t)) <= 0)
		return std::abs(cross(d1, d2, t)) / (d1 - d2).mag();
	return std::min((d1 - t).mag(), (d2 - t).mag());
}
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
//ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
ld area(const Polygon& H) {
	ld A = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) A += H[i] / H[(i + 1) % sz];
	return A * .5;
}
void norm(Polygon& H, const int& d = 1) {
	ld A = area(H);
	if (d == 1 && A < 0) std::reverse(H.begin(), H.end());
	else if (d == -1 && A > 0) std::reverse(H.begin(), H.end());
	return;
}
int inner_check(Pos H[], const int& sz, const Pos& p) {//concave
	int cnt = 0;
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
int inner_check_bi_search(Pos H[], const int& sz, const Pos& p) {//convex
	if (!sz) return -1;
	if (sz == 1) return p == H[0] ? 0 : -1;
	if (sz == 2) return on_seg_strong(H[0], H[1], p) ? 0 : -1;
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return -1;
	if (on_seg_strong(H[0], H[1], p) || on_seg_strong(H[0], H[sz - 1], p)) return 0;
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	if (cross(H[s], H[e], p) > 0) return 1;
	else if (on_seg_strong(H[s], H[e], p)) return 0;
	else return -1;
}
int inner_check_bi_search(const Polygon& H, const Pos& p) {//convex
	int sz = H.size();
	if (!sz) return -1;
	if (sz == 1) return p == H[0] ? 0 : -1;
	if (sz == 2) return on_seg_strong(H[0], H[1], p) ? 0 : -1;
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return -1;
	if (on_seg_strong(H[0], H[1], p) || on_seg_strong(H[0], H[sz - 1], p)) return 0;
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	if (cross(H[s], H[e], p) > 0) return 1;
	else if (on_seg_strong(H[s], H[e], p)) return 0;
	else return -1;
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
Polygon graham_scan(Polygon& C) {
	Polygon H;
	if (C.size() < 3) {
		std::sort(C.begin(), C.end());
		return C;
	}
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pos& p, const Pos& q) -> bool {
		int ret = ccw(C[0], p, q);
		if (!ret) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i]) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	return H;
}
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	std::sort(C.begin(), C.end());
	if (C.size() <= 2) { for (const Pos& p : C) H.push_back(p); }
	else {
		for (int i = 0; i < C.size(); i++) {
			while (H.size() > 1 && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
		int s = H.size() + 1;
		for (int i = C.size() - 1; i >= 0; i--) {
			while (H.size() > s && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
	}
	return H;
}
ld rotating_calipers(const Polygon& H) {
	int sz = H.size();
	if (sz < 2) return 0;
	if (sz == 2) return (H[0] - H[1]).mag();
	ld d = -1;
	for (int i = 0, j = 1; i < sz; i++) {
		while (ccw(H[i], H[(i + 1) % sz], H[j], H[(j + 1) % sz]) >= 0)
			j = (j + 1) % sz, d = std::max(d, (H[i] - H[j]).mag());
		d = std::max(d, (H[i] - H[j]).mag());
	}
	return d;
}
Polygon polygon_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
	Polygon qs;
	int n = ps.size();
	for (int i = 0; i < n; i++) {
		Pos p1 = ps[i], p2 = ps[(i + 1) % n];
		int d1 = ccw(b1, b2, p1), d2 = ccw(b1, b2, p2);
		if (d1 >= 0) qs.push_back(p1);
		if (d1* d2 < 0) qs.push_back(intersection(p1, p2, b1, b2));
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
	return ret;
}
struct Seg {
	Pos s, e, dir;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) { dir = e - s; }
	//bool operator < (const Seg& l) const { return s == l.s ? e < l.e : s < l.s; }
	bool inner(const Pos& p) const { return sign(dir / (p - s)) > 0; }
	friend bool parallel(const Seg& l0, const Seg& l1) { return zero(l0.dir / l1.dir); }
	friend bool same_dir(const Seg& l0, const Seg& l1) { return parallel(l0, l1) && l0.dir * l1.dir > 0; }
	friend Pos intersection_(const Seg& s1, const Seg& s2) {
		const Pos& p1 = s1.s, & p2 = s1.e;
		const Pos& q1 = s2.s, & q2 = s2.e;
		ld a1 = cross(q1, q2, p1);
		ld a2 = -cross(q1, q2, p2);
		return (p1 * a2 + p2 * a1) / (a1 + a2);
	}
	bool operator < (const Seg& l) const {
		if (same_dir(*this, l)) return l.inner(s);
		bool f0 = O < dir;
		bool f1 = O < l.dir;
		if (f0 != f1) return f1;
		return sign(dir / l.dir) > 0;
	}
	//bool operator == (const Seg& l) const { return s == l.s && e == l.e; }
	Seg operator + (const ld& d) const { Pos v = ~dir.unit(); return Seg(s - v * d, e - v * d); }
	Seg operator - (const ld& d) const { Pos v = ~dir.unit(); return Seg(s + v * d, e + v * d); }
	Seg operator += (const ld& d) { Pos v = ~dir.unit(); s -= v * d; e -= v * d; return *this; }
	Seg operator -= (const ld& d) { Pos v = ~dir.unit(); s += v * d; e += v * d; return *this; }
	Seg operator + (const Pos& v) const { return Seg(s + v, e + v); }
	Seg operator - (const Pos& v) const { return Seg(s - v, e - v); }
	Seg operator += (const Pos& v) { s += v; e += v; return *this; }
	Seg operator -= (const Pos& v) { s -= v; e -= v; return *this; }
	Seg operator * (const ld& d) const { return Seg(s, s + dir * d); }
	Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
typedef std::vector<Seg> Segs;
ld dot(const Seg& p, const Seg& q) { return dot(p.s, p.e, q.s, q.e); }
ld intersection(const Seg& s1, const Seg& s2, const int& f = 0) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == 2) return a1;
	if (f == 1) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
Segs half_plane_intersection(Segs& HP, const bool& srt = 1) {
	auto check = [&](Seg& u, Seg& v, Seg& w) -> bool {
		return w.inner(intersection_(u, v));
		};
	if (srt) std::sort(HP.begin(), HP.end());
	std::deque<Seg> dq;
	int sz = HP.size();
	for (int i = 0; i < sz; ++i) {
		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
		dq.push_back(HP[i]);
	}
	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
	sz = dq.size();
	if (sz < 3) return {};
	std::vector<Seg> HPI;
	for (int i = 0; i < sz; ++i) HPI.push_back(dq[i]);
	return HPI;
}
Segs half_plane_intersection(const Segs& P, const Segs& Q) {
	Segs HP(P.size() + Q.size());
	std::merge(P.begin(), P.end(), Q.begin(), Q.end(), HP.begin());
	return half_plane_intersection(HP, 0);
}
bool half_plane_intersection(Segs& HP, Segs& SHPI, Polygon& PHPI, const int& F = POS, const bool& srt = 1) {
	auto check = [&](Seg& u, Seg& v, Seg& w) -> bool {
		return w.inner(intersection_(u, v));
		};
	if (srt) std::sort(HP.begin(), HP.end());
	std::deque<Seg> dq;
	int sz = HP.size();
	for (int i = 0; i < sz; ++i) {
		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
		dq.push_back(HP[i]);
	}
	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
	sz = dq.size();
	if (sz < 3) return 0;
	if (F == POS) for (int i = 0; i < sz; ++i) PHPI.push_back(intersection_(dq[i], dq[(i + 1) % sz]));
	else if (F == SEG) for (int i = 0; i < sz; ++i) SHPI.push_back(dq[i]);
	return 1;
}
bool connectable(const Pos& u, const Pos& v, const Polygon& H, const bool& f = 1) {
	int sz = H.size();
	Pos m = (u + v) * .5;
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = H[i], p1 = H[(i + 1) % sz];
		if (intersect(p0, p1, u, v, WEAK)) return 0;
	}
	if (f == 1 && inner_check(H, m) == 2) return 0;
	if (f == 0 && inner_check(H, m) == 0) return 0;
	return 1;
}
ld C[LEN]; int vp;
struct Info {
	int i, t;
	ld c;
	Info(int i_ = 0, ld c_ = 0, int t_ = 0) : i(i_), c(c_), t(t_) {}
	bool operator < (const Info& x) const { return zero(c - x.c) ? i < x.i : c > x.c; }
};
std::vector<Info> G[LEN];
std::vector<Info> GW[LEN], GE[LEN];
int szw, sze;
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
ld dijkstra(std::vector<Info> G[], const int& v, const int& g) {
	for (int i = 0; i < LEN; i++) C[i] = INF;
	std::priority_queue<Info> Q; Q.push(Info(v, 0));
	C[v] = 0;
	while (Q.size()) {
		Info p = Q.top(); Q.pop();
		if (p.c > C[p.i]) continue;
		for (const Info& w : G[p.i]) {
			ld cost = p.c + w.c;
			if (C[w.i] > cost) {
				C[w.i] = cost;
				Q.push(Info(w.i, cost));
			}
		}
	}
	return C[g];
}
struct Event {
	Seg w, e;
	Event(Seg w_ = Seg(), Seg e_ = Seg()) : w(w_), e(e_) {
		if (w.s != w.e && e.s != e.e) {
			ld lo = 0, hi = 1;
			Pos v = ~w.dir;
			ld wlo = intersection(w, Seg(e.e, e.e + v), 2);
			ld whi = intersection(w, Seg(e.s, e.s + v), 2);
			ld elo = intersection(e, Seg(w.e, w.e + v), 2);
			ld ehi = intersection(e, Seg(w.s, w.s + v), 2);
			w = Seg(w.p(std::max(wlo, lo)), w.p(std::min(whi, hi)));
			e = Seg(e.p(std::min(ehi, hi)), e.p(std::max(elo, lo)));
		}
	}
};
ld dijkstra(
	const Pos& s, const Pos& t,
	const Event& we, const ld& m,
	const Polygon& W, const Polygon& E, const Polygon& H) {
	for (int i = 0; i < szw; i++) while (GW[i].back().t) GW[i].pop_back();
	for (int i = 0; i < sze; i++) while (GE[i].back().t) GE[i].pop_back();
	Pos w = we.w.p(m), e = we.e.p(m);
	int sz;
	if (connectable(w, s, H)) GW[0].push_back(Info(1, (s - w).mag()));
	sz = W.size();
	for (int i = 0; i < sz; i++)
		if (connectable(w, W[i], H)) {
			ld d = (w - W[i]).mag();
			GW[1].push_back(Info(i + 2, d, 1));
			GW[i + 2].push_back(Info(1, d, 1));
		}
	if (connectable(e, t, H)) GE[0].push_back(Info(1, (t - e).mag()));
	sz = E.size();
	for (int i = 0; i < sz; i++)
		if (connectable(e, E[i], H)) {
			ld d = (e - E[i]).mag();
			GE[1].push_back(Info(i + 2, d, 1));
			GE[i + 2].push_back(Info(1, d, 1));
		}
	ld dw = dijkstra(GW, 0, 1);
	ld de = dijkstra(GE, 0, 1);
	return dw + de;
}
ld ternary_search(
	const Pos& s, const Pos& t,
	const Event& we,
	const Polygon& W, const Polygon& E, const Polygon& H) {
	ld lo = 0, hi = 1;
	int c = 100;
	while (c--) {
		ld m1 = (lo + lo + hi) / 3;
		ld m2 = (lo + hi + hi) / 3;
		ld d1 = dijkstra(s, t, we, m1, W, E, H);
		ld d2 = dijkstra(s, t, we, m2, W, E, H);
		if (d1 < d2) lo = m1;
		else d2 = m2;
	}
	return (lo + hi) * .5;
}
#define BOJ
#ifdef BOJ
bool add_node() {

}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	//freopen("../../tests/aed/input/input0.txt", "r", stdin);
	//freopen("../../tests/aed/input/ret.txt", "w", stdout);
	Pos s, t; std::cin >> s >> t;
	std::cin >> N; Polygon W(N); for (Pos& p : W) std::cin >> p;
	if (!eq(W.back().y, 0)) std::reverse(W.begin(), W.end());
	std::cin >> M; Polygon E(M); for (Pos& p : E) std::cin >> p;
	if (!eq(E.back().y, 1000)) std::reverse(E.begin(), E.end());
	K = N + M; Polygon H = W; for (const Pos& p : E) H.push_back(p);
	//std::reverse(E.begin(), E.end());
	ld b = INF;
	Segs B;

	//connect land
	for (int i = 0; i < N; i++) {
		if (connectable(s, W[i], H)) {
			ld d = (s - W[i]).mag();
			GW[0].push_back(Info(i + 2, d, 0));
			GW[i + 2].push_back(Info(0, d, 0));
			G[0].push_back(Info(i + 2, d, 0));
			G[i + 2].push_back(Info(0, d, 0));
		}
		for (int j = i + 1; j < N; j++) {
			if (connectable(W[i], W[j], H)) {
				ld d = (W[i] - W[j]).mag();
				GW[i + 2].push_back(Info(j + 2, d, 0));
				GW[j + 2].push_back(Info(i + 2, d, 0));
				G[i + 2].push_back(Info(j + 2, d, 0));
				G[j + 2].push_back(Info(i + 2, d, 0));
			}
		}
	}
	for (int i = 0; i < M; i++) {
		if (connectable(t, E[i], H)) {
			ld d = (t - E[i]).mag();
			GE[0].push_back(Info(i + 2, d, 0));
			GE[i + 2].push_back(Info(0, d, 0));
			GE[1].push_back(Info(N + i + 4, d, 0));
			GE[N + i + 4].push_back(Info(1, d, 0));
		}
		for (int j = i + 1; j < M; j++) {
			if (connectable(E[i], E[j], H)) {
				ld d = (E[i] - E[j]).mag();
				GE[N + i + 4].push_back(Info(N + j + 4, d, 0));
				GE[N + j + 4].push_back(Info(N + i + 4, d, 0));
			}
		}
	}

	//connect bridge and get events
	for (int n = 0; n < N; n++) {
		const Pos& w0 = W[n], w1 = W[(n + 1) % N];
		Seg w01 = Seg(w0, w1);
		for (int m = 0; m < M; m++) {
			const Pos& e0 = E[m], e1 = E[(m + 1) % M];
			Seg e01 = Seg(e0, e1);
			ld d;
			d = (w0 - e0).mag();
			if (connectable(w0, e0, H, 0)) {
				if (eq(d, b)) B.push_back(Seg(w0, e0));
				else if (d < b) {
					B.clear();
					B.push_back(Seg(w0, e0));
				}
			}
			if (m < M - 1 && between(w0, w1, e0)) {
				d = cross(w0, w1, e0) / (w0 - w1).mag();
				if (sign(b - d) >= 0) {
					Pos v = ~(w1 - w0);
					ld x = intersection(w01, Seg(e0, e0 + v), 2);
					Pos p = w01.p(x);
					if (connectable(p, e0, H, 0)) {
						if (eq(d, b)) B.push_back(Seg(p, e0));
						else if (d < b) {
							B.clear();
							B.push_back(Seg(p, e0));
						}
					}
				}
			}
		}
	}
	return;
}
#else
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	int Q = 130;
	for (int q = 0; q <= Q; q++) {
		std::string Din = "../../tests/aed/input/input";
		std::string Dout = "../../tests/aed/output/output";
		std::string Qin = Din + std::to_string(q) + ".txt";
		std::string Qout = Dout + std::to_string(q) + ".txt";
		freopen(Qin.c_str(), "r", stdin);
		//freopen("../../tests/aed/input/ret.txt", "w", stdout);
		Pos u;
		std::cin >> N; Polygon H(N); std::cin >> u;
		for (Pos& p : H) std::cin >> p; norm(H);
		Vint I, J;
		for (int i = 0, i1; i < N; i++) {
			i1 = (i + 1) % N;
			Pos p0 = H[i], p1 = H[i1];
			bool vis;
			if (!ccw(u, p0, p1)) {
				if (dot(p0, p1, u) < 0) std::swap(p0, p1);
				vis = 1;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					if (intersect(q0, q1, p1, u)) { vis = 0; break; }
				}
			}
			else {
				if (ccw(u, p0, p1) < 0) std::swap(p0, p1);
				Seg s1 = Seg(p0, p1);
				Vrange V;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					Seg s2 = Seg(q0, q1);
					Range r = range(s1, s2, u);
					if (r.s.den != -1) V.push_back(r);
				}
				vis = 0;
				std::sort(V.begin(), V.end());
				V.push_back(Range(Frac(1), Frac(1)));
				int sz = V.size();
				Frac hi = Frac(0);
				for (const Range& r : V) {
					if (hi < r.s) { vis = 1; break; }
					else hi = std::max(hi, r.e);
				}
			}
			if (vis) I.push_back(i + 1);
		}
		//std::cout << I.size() << "\n";
		//for (int i : I) std::cout << i << " ";
		freopen(Qout.c_str(), "r", stdin);
		std::cin >> M;
		J.resize(M); for (int& j : J) std::cin >> j;
		if (I.size() != M) { std::cout << "fuck::\n"; continue; }
		bool f = 1;
		for (int j = 0; j < M; j++) {
			if (I[j] != J[j]) { f = 0; std::cout << "fuck::\n"; break; }
		}
		if (f) std::cout << "good::\n";
		else std::cout << "fuck::\n";
	}
	return;
}
#endif
int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068

/*

12
66207546 277078203
66231235 277049680
66236583 277045669
66270598 277030914
66302758 277017619
66330166 277008238
66376807 277002726
66421946 277004478
66139790 277216095
66150747 277172271
66169098 277129039
66185778 277105354

*/