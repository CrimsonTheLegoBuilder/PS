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
const int LEN = 25;
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

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

#define ABS 0
#define REL 1

//freopen("../../tests/", "r", stdin);
//freopen("../../tests/", "w", stdout);

int N, M, K, T, Q;
struct Pos {
	ld x, y;
	int i, j;
	Pos(ld x_ = 0, ld y_ = 0, int i_ = -1, int j_ = -1) : x(x_), y(y_), i(i_), j(j_) {  }
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	//bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	bool operator < (const Pos& p) const {//sort ccw
		bool f1 = zero(x) ? 0 < y : 0 < x;
		bool f2 = zero(p.x) ? 0 < p.y : 0 < p.x;
		if (f1 != f2) return f1;
		ld tq = *this / p;
		return tq > 0;
	}
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
	ld rad() const { return atan2l(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} V[LEN * LEN * 4]; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
std::istream& operator >> (std::istream& is, Polygon& P) { for (Pos& p : P) is >> p.x >> p.y; return is; }
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
//int inner_check(Pos H[], const int& sz, const Pos& p) {//concave
//	int cnt = 0;
//	for (int i = 0; i < sz; i++) {
//		Pos cur = H[i], nxt = H[(i + 1) % sz];
//		if (on_seg_strong(cur, nxt, p)) return 1;
//		if (zero(cur.y - nxt.y)) continue;
//		if (nxt.y < cur.y) std::swap(cur, nxt);
//		if (nxt.y - TOL < p.y || cur.y > p.y) continue;
//		cnt += ccw(cur, nxt, p) > 0;
//	}
//	return (cnt & 1) * 2;
//}
int inner_check(const Polygon& H, const Pos& p) {//concave
	int cnt = 0, sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (zero(cur.y - nxt.y)) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y <= p.y || cur.y > p.y) continue;
		//if (nxt.y - TOL < p.y || cur.y > p.y) continue;
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
Polygon convex_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
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
		ret = convex_cut(ret, b1, b2);
	}
	return ret;
}
//struct Seg {
//	Pos s, e, dir;
//	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) { dir = e - s; }
//	//bool operator < (const Seg& l) const { return s == l.s ? e < l.e : s < l.s; }
//	bool inner(const Pos& p) const { return sign(dir / (p - s)) > 0; }
//	friend bool parallel(const Seg& l0, const Seg& l1) { return zero(l0.dir / l1.dir); }
//	friend bool same_dir(const Seg& l0, const Seg& l1) { return parallel(l0, l1) && l0.dir * l1.dir > 0; }
//	friend Pos intersection_(const Seg& s1, const Seg& s2) {
//		const Pos& p1 = s1.s, & p2 = s1.e;
//		const Pos& q1 = s2.s, & q2 = s2.e;
//		ld a1 = cross(q1, q2, p1);
//		ld a2 = -cross(q1, q2, p2);
//		return (p1 * a2 + p2 * a1) / (a1 + a2);
//	}
//	bool operator < (const Seg& l) const {
//		if (same_dir(*this, l)) return l.inner(s);
//		bool f0 = O < dir;
//		bool f1 = O < l.dir;
//		if (f0 != f1) return f1;
//		return sign(dir / l.dir) > 0;
//	}
//	//bool operator == (const Seg& l) const { return s == l.s && e == l.e; }
//	Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
//	ld green(const ld& lo = 0, const ld& hi = 1) const {
//		ld d = hi - lo;
//		ld ratio = (lo + hi) * .5;
//		Pos m = p(ratio);
//		return m.y * d * (s.x - e.x);
//	}
//};
//typedef std::vector<Seg> Segs;
//ld dot(const Seg& p, const Seg& q) { return dot(p.s, p.e, q.s, q.e); }
//ld intersection(const Seg& s1, const Seg& s2, const bool& f = 0) {
//	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
//	ld det = (q2 - q1) / (p2 - p1);
//	if (zero(det)) return -1;
//	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
//	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
//	if (f == 1) return fit(a1, 0, 1);
//	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
//	return -1;
//}
//Segs half_plane_intersection(Segs& HP, const bool& srt = 1) {
//	auto check = [&](Seg& u, Seg& v, Seg& w) -> bool {
//		return w.inner(intersection_(u, v));
//		};
//	if (srt) std::sort(HP.begin(), HP.end());
//	std::deque<Seg> dq;
//	int sz = HP.size();
//	for (int i = 0; i < sz; ++i) {
//		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
//		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
//		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
//		dq.push_back(HP[i]);
//	}
//	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
//	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
//	sz = dq.size();
//	if (sz < 3) return {};
//	std::vector<Seg> HPI;
//	for (int i = 0; i < sz; ++i) HPI.push_back(dq[i]);
//	return HPI;
//}
//Segs half_plane_intersection(const Segs& P, const Segs& Q) {
//	Segs HP(P.size() + Q.size());
//	std::merge(P.begin(), P.end(), Q.begin(), Q.end(), HP.begin());
//	return half_plane_intersection(HP, 0);
//}
//struct Circle {
//	Pos c;
//	int r;
//	Circle(Pos c_ = Pos(), int r_ = 0) : c(c_), r(r_) {}
//	bool operator == (const Circle& q) const { return c == q.c && r == q.r; }
//	bool operator != (const Circle& q) const { return !(*this == q); }
//	bool operator < (const Circle& q) const { return c == q.c ? r < q.r : c < q.c; }
//	//bool operator < (const Circle& q) const { return r < q.r && (c - q.c).mag() + r < q.r + TOL; }
//	bool outside(const Circle& q) const { return sign((c - q.c).Euc() - sq((ll)r + q.r)) >= 0; }
//	Circle operator + (const Circle& q) const { return { c + q.c, r + q.r }; }
//	Circle operator - (const Circle& q) const { return { c - q.c, r - q.r }; }
//	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
//	ld rad(const Pos& p) const { return (p - c).rad(); }
//	ld area(const ld& lo, const ld& hi) const { return (hi - lo) * r * r * .5; }
//	ld green(const ld& lo, const ld& hi) const {
//		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
//		ld fan = area(lo, hi);
//		Pos m = c + (s + e) * r * (ld).5;
//		ld tz = (cos(lo) - cos(hi)) * m.y * r;
//		return fan + tz - (s / e) * r * r * (ld).5;
//	}
//	ld H(const ld& th) const { return sin(th) * c.x + cos(th) * c.y + r; }//coord trans | check right
//	//bool operator < (const Pos& p) const { return r < (c - p).mag(); }
//	bool operator < (const Pos& p) const { return sign(r - (c - p).mag()) < 0; }
//	bool operator > (const Pos& p) const { return r > (c - p).mag(); }
//	bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
//	friend std::istream& operator >> (std::istream& is, Circle& c) { is >> c.c >> c.r; return is; }
//	friend std::ostream& operator << (std::ostream& os, const Circle& c) { os << c.c << " " << c.r; return os; }
//} INVAL = { { 0, 0 }, -1 };
//bool cmpr(const Circle& p, const Circle& q) { return p.r > q.r; }//sort descending order
//Vld intersections(const Circle& a, const Circle& b) {
//	Pos ca = a.c, cb = b.c;
//	Pos vec = cb - ca;
//	ll ra = a.r, rb = b.r;
//	ld distance = vec.mag();
//	ld rd = vec.rad();
//	if (vec.Euc() > sq(ra + rb) + TOL) return {};
//	if (vec.Euc() < sq(ra - rb) - TOL) return {};
//	ld X = (ra * ra - rb * rb + vec.Euc()) / (2 * distance * ra);
//	if (X < -1) X = -1;
//	if (X > 1) X = 1;
//	ld h = acos(X);
//	Vld ret = {};
//	ret.push_back(norm(rd + h));
//	if (zero(h)) return ret;
//	ret.push_back(norm(rd - h));
//	return ret;
//}
//Vld circle_line_intersections(const Seg& l, const Circle& q, const int& t = LINE) {
//	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
//	Pos s = l.s, e = l.e;
//	Pos vec = e - s;
//	Pos OM = s - q.c;
//	ld a = vec.Euc();
//	ld b = vec * OM;
//	ld c = OM.Euc() - q.r * q.r;
//	ld J = b * b - a * c;
//	if (J < -TOL) return {};
//	ld det = sqrt(std::max((ld)0, J));
//	ld lo = (-b - det) / a;
//	ld hi = (-b + det) / a;
//	Vld ret;
//	if (t == LINE) {
//		if (0 < lo && lo < 1) ret.push_back(lo);
//		if (zero(det)) return ret;
//		if (0 < hi && hi < 1) ret.push_back(hi);
//	}
//	else {//circle
//		auto the = [&](ld rt) { return q.rad(s + (e - s) * rt); };
//		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
//		if (zero(det)) return ret;
//		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
//	}
//	return ret;
//}
//Vld circle_line_intersections(const Pos& s, const Pos& e, const Circle& q, const bool& f = 0) {
//	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
//	Pos vec = e - s;
//	Pos OM = s - q.c;
//	ld a = vec.Euc();
//	ld b = vec * OM;
//	ld c = OM.Euc() - q.r * q.r;
//	ld J = b * b - a * c;
//	if (J < -TOL) return {};
//	ld det = sqrt(std::max((ld)0, J));
//	ld lo = (-b - det) / a;
//	ld hi = (-b + det) / a;
//	Vld ret;
//	if (f) {
//		if (0 < hi && hi < 1) ret.push_back(hi);
//		if (zero(det)) return ret;
//		if (0 < lo && lo < 1) ret.push_back(lo);
//	}
//	else {
//		auto the = [&](ld rt) { return q.rad(s + (e - s) * rt); };
//		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
//		if (zero(det)) return ret;
//		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
//	}
//	return ret;
//}
Polygon circle_line_intersection(const Pos& o, const ld& r, const Pos& p1, const Pos& p2) {
	ld d = dist(p1, p2, o);
	if (std::abs(d) > r) return {};
	Pos vec = p2 - p1;
	Pos m = intersection(p1, p2, o, o + ~vec);
	ld distance = vec.mag();
	ld ratio = sqrt(r * r - d * d);
	Pos m1 = m - vec * ratio / distance;
	Pos m2 = m + vec * ratio / distance;
	if (dot(p1, p2, m1, m2) < 0) std::swap(m1, m2);
	return { m1, m2 };//p1->p2
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
//ld circle_cut(const Circle& c, const Pos& p1, const Pos& p2) {
//	Pos v1 = p1 - c.c, v2 = p2 - c.c;
//	ld r = c.r;
//	std::vector<Pos> inx = circle_line_intersection(O, r, v1, v2);
//	if (inx.empty()) return r * r * rad(v1, v2) * .5;
//	Pos m1 = inx[0], m2 = inx[1];
//	bool d1 = dot(m1, v1, m2) > -TOL, d2 = dot(m1, v2, m2) > -TOL;
//	if (d1 && d2) return (v1 / v2) * .5;
//	else if (d1) return (v1 / m2 + r * r * rad(m2, v2)) * .5;
//	else if (d2) return (r * r * rad(v1, m1) + m1 / v2) * .5;
//	else if (dot(v1, m1, v2) > 0 && dot(v1, m2, v2) > 0)
//		return (r * r * (rad(v1, m1) + rad(m2, v2)) + m1 / m2) * .5;
//	else return (r * r * rad(v1, v2)) * .5;
//}
//struct Arc {
//	ld lo, hi;
//	Arc(ld l_ = 0, ld h_ = 0) : lo(l_), hi(h_) {}
//	bool operator < (const Arc& a) const { return zero(lo - a.lo) ? hi < a.hi : lo < a.lo; }
//	inline friend std::istream& operator >> (std::istream& is, Arc& a) { is >> a.lo >> a.hi; return is; }
//	inline friend std::ostream& operator << (std::ostream& os, const Arc& a) { os << a.lo << " " << a.hi; return os; }
//};
ld C[LEN * LEN * 4]; int vp;
struct Info {
	int i;
	ld c;
	Info(int i_ = 0, ld c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return c > x.c; }
};
std::vector<Info> G[LEN * LEN * 4];
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
Polygon RV[LEN];//revolve
ld get_theta(const Pos& d1, const Pos& d2, const ld& r) { return asin(r / (d1 - d2).mag()); }
Pos mid(const Pos& d1, const Pos& d2) { return (d1 + d2) * .5; }
bool between(const Pos& d1, const Pos& d2, const Pos& q) { return sign(dot(d1, d2, q)) <= 0 && sign(dot(d2, d1, q)) <= 0; }
bool close(const Pos& d, const Pos& q, const ld& r) { return (d - q).mag() < r; }
bool close(const Pos& d1, const Pos& d2, const Pos& q, const ld& r) { dist(d1, d2, q, ABS) < r; }
Pos rotate(const Pos& p, const Pos& pv, const ld& t, const int& i) {
	Pos v = p - pv;
	ld rat = cos(t);
	Pos q = v.rot(t) * rat; q += pv; q.i = i;
	return q;
}
bool circle_is_ok(const Pos& c, const ld& r, const Polygon& P) {
	int sz = P.size();
	//std::cout << "c:: " << c << "\n";
	//std::cout << "r:: " << r << "\n";
	for (int i = 0; i < sz; i++) {
		//std::cout << "P[" << i << "]:: " << P[i] << ", P[" << i + 1 << "]:: " << P[(i + 1) % sz] << "\n";
		//std::cout << "dist:: " << dist(P[i], P[(i + 1) % sz], c, ABS) << "\n";
		if (sign(dist(P[i], P[(i + 1) % sz], c, ABS) - r + TOL) < 0) return 0;
	}
	if (inner_check(P, c)) return 0;
	return 1;
}
bool connectable(const Pos& s, const Pos& e, const ld& r, const Polygon& P) {
	if (s == e) return 1;
	Pos v = ~(e - s).unit() * r;
	Polygon B = { s + v, s - v, e - v, e + v };
	ld a = area(sutherland_hodgman(P, B));
	return zero(a);
}
void connect_node(const int& n1, const int& n2, const Polygon& P, const ld& r) {
	Pos d1 = V[n1], d2 = V[n2];
	if (d1.i != d2.i && connectable(d1, d2, r, P)) {
		G[n1].push_back({ n2, (d1 - d2).mag() });
		G[n2].push_back({ n1, (d1 - d2).mag() });
	}
	return;
}
void connect_seg(const Polygon& P, const ld& r) {
	for (int i = 0; i < vp; i++)
		for (int j = i + 1; j < vp; j++)
			connect_node(i, j, P, r);
	return;
}
void connect_arc(const Polygon& P, const ld& r) {
	int psz = P.size();
	for (int i = 0; i < psz; i++) {
		int sz = RV[i].size();
		for (int j = 0; j < sz; j++) RV[i][j] -= P[i];
		std::sort(RV[i].begin(), RV[i].end());
		for (int j = 0; j < sz; j++) RV[i][j] += P[i];
		for (int j = 0; j < sz; j++) {
			Pos cur = RV[i][j], nxt = RV[i][(j + 1) % sz];
			if (cur == nxt) {
				G[cur.j].push_back({ nxt.j, 0 });
				G[nxt.j].push_back({ cur.j, 0 });
				continue;
			}
			bool f0 = 1;
			for (int k = 0; k < psz; k++) {
				if (i != k) {
					bool f1 = inside(nxt, P[i], cur, P[k]);
					bool f2 = sign((P[i] - P[k]).mag() - r * 2) < 0;
					bool f3 = sign((cur - P[k]).mag() - r) < 0 
						|| sign((nxt - P[k]).mag() - r) < 0;
					if ((f1 && f2) || f3) { f0 = 0; break; }
				}
				const Pos& p0 = P[(k - 1 + psz) % psz], & p1 = P[k];
				Polygon inx = circle_seg_intersection(P[i], r * 2, p0, p1);
				for (const Pos& p : inx) {
					if (inside(nxt, P[i], cur, p)) {
						ld d = dist(p0, p1, P[i]);
						if (d < r * 2) { f0 = 0; break; }
					}
				}
			}
			if (f0) {
				Pos p = P[i];
				ld t1 = norm((cur - p).rad());
				ld t2 = norm((nxt - p).rad());
				ld t = norm(t2 - t1);
				ld rd = r * t;
				G[cur.j].push_back({ nxt.j, rd });
				G[nxt.j].push_back({ cur.j, rd });
			}
		}
	}
	return;
}
void pos_init(const Pos& s, const Pos& e, const Polygon& P, const ld& r) {
	vp = 0;
	V[vp++] = s;
	V[vp++] = e;
	int sz = P.size();
	for (int i = 0; i < sz; i++) {//tangent from S || E
		ld t1 = get_theta(s, P[i], r);
		ld t2 = get_theta(e, P[i], r);
		V[vp++] = rotate(P[i], s, t1, i);
		if (!circle_is_ok(V[vp - 1], r, P)) vp--;
		V[vp++] = rotate(P[i], s, -t1, i);
		if (!circle_is_ok(V[vp - 1], r, P)) vp--;
		V[vp++] = rotate(P[i], e, t2, i);
		if (!circle_is_ok(V[vp - 1], r, P)) vp--;
		V[vp++] = rotate(P[i], e, -t2, i);
		if (!circle_is_ok(V[vp - 1], r, P)) vp--;
	}
	for (int i = 0; i < sz; i++) {
		for (int j = i + 1; j < sz; j++) {
			if (i == j) continue;
			Pos v = ~(V[j] - V[i]).unit();
			V[vp++] = P[i] + v * r;
			if (!circle_is_ok(V[vp - 1], r, P)) vp--;
			V[vp++] = P[j] + v * r;
			if (!circle_is_ok(V[vp - 1], r, P)) vp--;
			V[vp++] = P[i] - v * r;
			if (!circle_is_ok(V[vp - 1], r, P)) vp--;
			V[vp++] = P[j] - v * r;
			if (!circle_is_ok(V[vp - 1], r, P)) vp--;
			ld d = (V[j] - V[i]).mag();
			if (d > r * 2) {//tangent from m
				Pos m = mid(P[i], P[j]);
				ld t = get_theta(m, P[i], r);
				V[vp++] = rotate(P[i], m, t, i);
				if (!circle_is_ok(V[vp - 1], r, P)) vp--;
				V[vp++] = rotate(P[j], m, t + PI, j);
				if (!circle_is_ok(V[vp - 1], r, P)) vp--;
				V[vp++] = rotate(P[i], m, -t, i);
				if (!circle_is_ok(V[vp - 1], r, P)) vp--;
				V[vp++] = rotate(P[j], m, -t - PI, j);
				if (!circle_is_ok(V[vp - 1], r, P)) vp--;
			}
		}
	}
	std::cout << "vp:: " << vp << "\n";
	assert(vp < 2500);
	for (int i = 0; i < vp; i++) G[i].clear();
	for (int i = 0; i < 20; i++) RV[i].clear();
	for (int i = 2; i < vp; i++) { V[i].j = i; RV[V[i].i].push_back(V[i]); }
	std::cout << "DEBUG::\n";
	for (int i = 0; i < vp; i++) {
		std::cout << "V[" << i << "] = (";
		std::cout << V[i].x << ", " << V[i].y << ")\n";
	}
	std::cout << "DEBUG::\n";
	return;
}
bool query() {
	std::cin >> N;
	if (!N) return 0;
	//Polygon P(N); for (Pos& p : P) std::cin >> p; norm(P);
	Polygon P(N); std::cin >> P; norm(P); for (int i = 0; i < N; i++) P[i].i = i;
	Pos s, e; std::cin >> s >> e; s.i = -1; e.i = -2;
	ld r = 1;
	pos_init(s, e, P, r);
	connect_seg(P, r);
	connect_arc(P, r);
	ld d = dijkstra(0, 1);
	assert(d < INF);
	std::cout << d << "\n";
	return 1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	while (query());
	return;
}
int main() { solve(); return 0; }//boj22801