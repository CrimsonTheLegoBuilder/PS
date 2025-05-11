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
const int LEN = 25;
inline int sign(const ld& x) { return x <= -TOL ? -1 : x >= TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ll sq(const ll& x) { return x * x; }
inline ld sq(const ld& x) { return x * x; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }

#define INSIDE 0
#define MEET 1
#define OUTSIDE 2
#define LO x
#define HI y
#define LINE 1
#define CIRCLE 2
#define STRONG 0
#define WEAK 1
#define ABS 0
#define REL 1

ld heron(const ld& a, const ld& b, const ld& c) {
	ld s = (a + b + c) / 2;
	ld A = sqrt(s * (s - a) * (s - b) * (s - c));
	return A;
}
ld heron(const Vld& v) { assert(3 == v.size()); return heron(v[0], v[1], v[2]); }
int T, N, M;
bool F[3];
Vld Q;
Vint sts;
struct Pos {
	ld x, y;
	int i, j;
	Pos(ld x_ = 0, ld y_ = 0, int i_ = -1, int j_ = -1) : x(x_), y(y_), i(i_), j(j_) {}
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
typedef std::set<Pos> MapPos;
std::istream& operator >> (std::istream& is, Polygon& P) { for (Pos& p : P) is >> p.x >> p.y; return is; }
bool cmpr(const Pos& p, const Pos& q) {
	bool f1 = O < p;
	bool f2 = O < q;
	if (f1 != f2) return f1;
	ld tq = p / q;
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
	return ret;
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
	Pos p(const ld& rt) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
ld intersection(const Seg& s1, const Seg& s2, const bool& f = STRONG) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == WEAK) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
struct Circle {
	Pos c;
	ll r;
	Circle(Pos c_ = Pos(), ll r_ = 0) : c(c_), r(r_) {}
	bool operator > (const Pos& p) const { return sign(r - (c - p).mag()) > 0; }
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld rad(const Pos& p) const { return (p - c).rad(); }
	ld area(const ld& lo, const ld& hi) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
	friend std::istream& operator >> (std::istream& is, Circle& p) { is >> p.c.x >> p.c.y >> p.r; return is; }
	friend std::ostream& operator << (std::ostream& os, const Circle& p) { os << p.c.x << " " << p.c.y << " " << p.r; return os; }
} C[3];
Vld intersections(const Circle& a, const Circle& b) {
	Pos ca = a.c, cb = b.c;
	Pos vec = cb - ca;
	ld ra = a.r, rb = b.r;
	ld distance = vec.mag();
	ld rd = vec.rad();
	if (vec.Euc() > sq(ra + rb) + TOL) return {};
	if (vec.Euc() < sq(ra - rb) - TOL) return {};
	ld X = (ra * ra - rb * rb + vec.Euc()) / (2 * distance * ra);
	if (X < -1) X = -1;
	if (X > 1) X = 1;
	ld h = acos(X);
	Vld ret = {};
	ret.push_back(norm(rd - h));
	if (zero(h)) return ret;
	ret.push_back(norm(rd + h));
	return ret;
}
struct Sphere {
	ll x, y, z, r;
	Sphere(ll x_ = 0, ll y_ = 0, ll z_ = 0, ll r_ = 0) : x(x_), y(y_), z(z_), r(r_) {}
	bool operator < (const Sphere& q) const { return r > q.r; }
	Sphere operator - (const Sphere& q) const { return { x - q.x, y - q.y, z - q.z, 0 }; }
	ll operator * (const Sphere& q) const { return x * q.x + y * q.y + z * q.z; }
	ld vol() const { return (4. / 3) * PI * r * r * r; }
	ld vol(const ld& h) const { return PI * h * h * (3 * r - h) / 3; }
	ld surf() const { return PI * 4 * r * r; }
	ld surf(const ld& h) const { return PI * 2 * r * h; }
} S[3], SP[LEN * 3];
ll Euc(const Sphere& p, const Sphere& q) { return sq(p.x - q.x) + sq(p.y - q.y) + sq(p.z - q.z); }
ld mag(const Sphere& p, const Sphere& q) { return sqrtl(Euc(p, q)); }
ld rad(const Sphere& a, const Sphere& b, const Sphere& c) {
	ld dab = mag(a, b);
	ld dac = mag(a, c);
	ll det = (b - a) * (c - a);
	ld proj = det / dab;
	ld ret = fit(proj / dac, -1, 1);
	return acos(ret);
}
int meet(const Sphere& p, const Sphere& q) {
	ll dist = Euc(p, q);
	ll rout = sq(p.r + q.r);
	if (dist >= rout) return OUTSIDE;
	ll rin = sq(p.r - q.r);
	if (dist <= rin) return INSIDE;
	return MEET;
}
ld cos_2nd(const ld& a, const ld& b, const ld& c) {
	ld num = sq(a) + sq(b) - sq(c);
	ld den = 2 * a * b;
	ld t = num / den;
	return std::abs(acosl(std::min(std::max(t, -(ld)1.0), (ld)1.0)));
}
void spherical_triangle_angles(const ld& a, const ld& b, const ld& c, ld& A_, ld& B_, ld& C_) {
	A_ = acos(fit((cos(a) - cos(b) * cos(c)) / (sin(b) * sin(c)), -1, 1));
	B_ = acos(fit((cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c)), -1, 1));
	C_ = acos(fit((cos(c) - cos(a) * cos(b)) / (sin(a) * sin(b)), -1, 1));
	return;
}
ld area(const ld& a, const ld& b, const ld& c, const ll& r, const ld& t) {
	ld A_, B_, C_;
	if (a >= PI) {
		spherical_triangle_angles(a * .5, t, c, A_, B_, C_);
		return r * r * (A_ + B_ + C_ - PI) * 2;
	}
	spherical_triangle_angles(a, b, c, A_, B_, C_);
	return r * r * (A_ + B_ + C_ - PI);
}
//ld two_union(const Sphere& a, const Sphere& b) {
//	int f = meet(a, b);
//	if (f == OUTSIDE) return a.vol() + b.vol();
//	if (f == INSIDE) { return a.r >= b.r ? a.vol() : b.vol(); }
//	ld d = mag(a, b);
//	Pos ca = Pos(0, 0);
//	Pos cb = Pos(d, 0);
//	Circle aa = Circle(ca, a.r);
//	Circle bb = Circle(cb, b.r);
//	Vld inxs = intersections(aa, bb);
//	assert(inxs.size());
//	Pos p = aa.p(inxs[0]);
//	ld ha = a.r - p.x;
//	ld hb = p.x - d + b.r;
//	return a.vol(a.r + a.r - ha) + b.vol(b.r + b.r - hb);
//}
ld cone_vol(const ld& r, const ld& h, const ld& h1, const ld& t = 2 * PI) {
	ld v = sq(r) * h * PI / 3;
	ld v0 = sq(r) * (h - h1) * PI / 3;
	return (v - v0) * (t / (2 * PI));
}
ld tri_pil_vol(const Vld& T, const Vld& H) {
	assert(T.size() == 3);
	assert(H.size() == 3);
	ld A = heron(T);
	return (H[0] + H[1] + H[2]) * A / 3;
}
ld two_union(const Sphere& a, const Sphere& b) {
	int f = meet(a, b);
	if (f == INSIDE) { return a.r >= b.r ? a.vol() : b.vol(); }
	ld d = mag(a, b);
	if (a.r == b.r) return a.vol() + sq(a.r) * PI * d;
	Pos ca = Pos(0, 0);
	Pos cb = Pos(d, 0);
	Circle aa = Circle(ca, a.r);
	Circle bb = Circle(cb, b.r);
	Pos v = aa.c - bb.c;
	ld hu = d;
	ld gu = a.r - b.r;
	ld go = sqrt(sq(hu) - sq(gu));
	ld t = norm(atan2(go, gu));
	ld ha = a.r * (1 + cos(t));
	ld hb = b.r * (1 - cos(t));
	ld r0 = a.r * std::abs(sin(t));
	ld r1 = b.r * std::abs(sin(t));
	if (r0 < r1) std::swap(r0, r1);
	ld t0 = std::abs(PI * .5 - t);
	ld h = r0 / tan(t0);
	ld h1 = go / cos(t0);
	return a.vol(a.r + a.r - ha) + b.vol(b.r + b.r - hb) + cone_vol(r0, h, h1);
}
Pos tangent(const Circle& c1, const Circle& c2) {
	Pos v = c2.c - c1.c;
	ld t0 = v.rad();
	Pos r1 = v.unit() * c1.r;
	Pos r2 = v.unit() * c2.r;
	ld hu = v.mag();
	ld gu = c1.r - c2.r;
	ld go = sqrt(sq(hu) - sq(gu));
	ld t = norm(atan2(go, gu));
	return Pos(t0 - t, t0 + t);
}
Polygon tangent_convex_hull(const Circle& c1, const Circle& c2) {
	Pos v = c2.c - c1.c;
	Pos r1 = v.unit() * c1.r;
	Pos r2 = v.unit() * c2.r;
	ld hu = v.mag();
	ld gu = c1.r - c2.r;
	ld go = sqrt(sq(hu) - sq(gu));
	ld t = norm(atan2(go, gu));
	Pos p1 = c1.c + r1.rot(t);
	Pos p2 = c1.c + r1.rot(-t);
	Pos p3 = c2.c + r2.rot(t);
	Pos p4 = c2.c + r2.rot(-t);
	Polygon C = { p1, p2, p3, p4 };
	return graham_scan(C);
}
bool inner_check(const Polygon& B, const Circle& c) {
	int sz = B.size();
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = B[i], & p1 = B[(i + 1) % sz];
		ld d = dist(p0, p1, c.c);
		if (d < c.r) return 0;
	}
	return 1;
}
bool intersect(const Circle& c, const Pos& hp1, const Pos& hp2) {
	Pos p1 = c.p(hp1.LO);
	Pos p2 = c.p(hp1.HI);
	Pos q1 = c.p(hp2.LO);
	Pos q2 = c.p(hp2.HI);
	return intersect(p1, p2, q1, q2);
}
struct Pos3D {
	ld x, y, z;
	Pos3D(ld x_ = 0, ld y_ = 0, ld z_ = 0) : x(x_), y(y_), z(z_) {}
	bool operator == (const Pos3D& p) const { return zero(x - p.x) && zero(y - p.y) && zero(z - p.z); }
	bool operator != (const Pos3D& p) const { return !zero(x - p.x) || !zero(y - p.y) || !zero(z - p.z); }
	bool operator < (const Pos3D& p) const { return zero(x - p.x) ? zero(y - p.y) ? z < p.z : y < p.y : x < p.x; }
	ld operator * (const Pos3D& p) const { return x * p.x + y * p.y + z * p.z; }
	Pos3D operator / (const Pos3D& p) const { return { y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x }; }
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	Pos3D operator * (const ld& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ld& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ld& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ld& n) { x /= n; y /= n; z /= n; return *this; }
	ld Euc() const { return x * x + y * y + z * z; }
	ld mag() const { return sqrtl(Euc()); }
	ld lon() const { return atan2(y, x); }
	ld lat() const { return atan2(z, sqrtl(x * x + y * y)); }
	Pos3D operator - () const { return { -x, -y, -z }; }
	Pos3D unit() const { return *this / mag(); }
	Pos3D norm(const Pos3D& p) const { return (*this / p).unit(); }
	Pos3D rodrigues_rotate(const ld& th, const Pos3D& axis) const {
		ld s = sin(th), c = cos(th); Pos3D n = axis.unit();
		return n * (*this * n) * (1 - c) + (*this * c) + (n / *this) * s;
	}
	friend ld rad(const Pos3D& p1, const Pos3D& p2) { return atan2l((p1 / p2).mag(), p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z; return os; }
};
typedef std::vector<Pos3D> Polyhedron;
const Pos3D O3D = { 0, 0, 0 };
const Pos3D X_axis = { 1, 0, 0 };
const Pos3D Y_axis = { 0, 1, 0 };
const Pos3D Z_axis = { 0, 0, 1 };
const Pos3D MAXP3D = { INF, INF, INF };
Pos3D pos3d(const Pos& p) { return Pos3D(p.x, p.y, 0); }
Pos3D pos3d(const Circle& c) { return Pos3D(c.c.x, c.c.y, c.r); }
Pos3D S2C(const ld& lon, const ld& lat) {//Spherical to Cartesian
	ld phi = lon * PI / 180;
	ld the = lat * PI / 180;
	return Pos3D(cos(phi) * cos(the), sin(phi) * cos(the), sin(the));
}
Pos3D point(const Pos3D Xaxis, const Pos3D Yaxis, const ld& th) { return Xaxis * cos(th) + Yaxis * sin(th); }
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) { return norm(atan2(Yaxis * p, Xaxis * p)); }
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
ld dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
ld area(const std::vector<Pos3D>& H, const Pos3D& norm) {
	if (H.size() < 3) return 0;
	ld ret = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos3D& cur = H[i], nxt = H[(i + 1) % sz];
		ret += cross(H[0], cur, nxt) * norm / norm.mag();
	}
	return std::abs(ret * .5);
}
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) > 0; }
struct Line3D {
	Pos3D dir, p0;
	Line3D(Pos3D dir_ = Pos3D(0, 0, 0), Pos3D p0_ = Pos3D(0, 0, 0)) : dir(dir_), p0(p0_) {}
};
Line3D L(const Pos3D& p1, const Pos3D& p2) { return { p2 - p1, p1 }; }
Line3D line(const Pos3D& p1, const Pos3D& p2) { return { p2 - p1, p1 }; }
struct Plane {
	ld a, b, c, d;
	Plane(ld a_ = 0, ld b_ = 0, ld c_ = 0, ld d_ = 0) : a(a_), b(b_), c(c_), d(d_) {}
	Pos3D norm() const { return Pos3D(a, b, c); };
	Plane operator + (const ld& n) const { return { a, b, c, d + n }; }
	Plane operator - (const ld& n) const { return { a, b, c, d - n }; }
	Plane& operator += (const ld& n) { d += n; return *this; }
	Plane& operator -= (const ld& n) { d -= n; return *this; }
	friend std::istream& operator >> (std::istream& is, Plane& f) { is >> f.a >> f.b >> f.c >> f.d; return is; }
	friend std::ostream& operator << (std::ostream& os, const Plane& f) { os << f.a << " " << f.b << " " << f.c << " " << f.d; return os; }
} knife;
Plane plane(const Pos3D& p, const ld& n) { return Plane(p.x, p.y, p.z, n); }
ld dist(const Plane& s, const Pos3D& p) { return (s.norm() * p + s.d) / s.norm().mag(); }
Pos3D intersection(const Plane& S, const Line3D& l) {
	ld det = S.norm() * l.dir;
	if (zero(det)) return { INF, INF, INF };
	ld t = -((S.norm() * l.p0 + S.d) / det);
	return l.p0 + (l.dir * t);
}
Pos3D intersection(const Plane& S, const Pos3D& p1, const Pos3D& p2, const bool& f = 0) {
	Line3D l = line(p1, p2);
	Pos3D inx = intersection(S, l);
	if (f && !on_seg_strong(p1, p2, inx)) return { INF, INF, INF };
	return inx;
}
struct Planar {
	Pos3D n, p0;
	Planar(Pos3D n_ = Pos3D(0, 0, 0), Pos3D p0_ = Pos3D(0, 0, 0)) : n(n_), p0(p0_) {}
};
Planar planar(const Pos3D& p1, const Pos3D& p2, const Pos3D& p3) {
	Pos3D norm = (p2 - p1) / (p3 - p2);
	return Planar(norm, p1);
}
ld rad(const Planar& p1, const Planar& p2) { return norm(PI - rad(p1.n, p2.n)); }
ld volume(const ll& r, const Polygon& hp) {
	int sz = hp.size();
	assert(sz);
	if (sz == 1) {
		ld t = norm(hp[0].HI - hp[0].LO) * .5;
		ld d = r * cosl(t);
		ld h = r - d;
		return Sphere(0, 0, 0, r).vol(r + r - h);
	}
	auto inside = [&](const Pos& p, const ld& t) -> bool {
		if (p.LO < p.HI) {
			if (p.LO < t && t < p.HI) return 1;
		}
		else {//(p.LO > p.HI)
			if (p.LO < t || t < p.HI) return 1;
		}
		return 0;
		};
	auto outer_check = [&](const Pos& p, const Pos& q) -> bool {
		return !inside(p, q.LO) && !inside(p, q.HI);
		};
	auto the = [&](const ld& dd, const ld& rr) -> ld {
		ld w = dd / rr;
		assert(std::abs(dd) <= std::abs(rr));
		return norm(acosl(fit(w, -1, 1))) * 2;
		};
	auto area_ = [&](const ld& rr, const ld& t) -> ld {
		ld fan = std::abs(rr * rr * (2 * PI - t) * .5);
		ld z = t * .5;
		ld tri = rr * sin(z) * rr * cos(z);
		tri = std::abs(tri);
		if (t < PI) return fan + tri;
		return fan - tri;
		};
	auto cone_vol = [&](const ld& rr, const ld& t, const ld& h) -> ld {
		ld rat = area_(rr, t) / (rr * rr * PI);
		ld vol = rr * rr / 3 * rat * PI * h;
		return vol;
		};
	Circle c = Circle(Pos(), r);
	assert(sz == 2);
	Pos u = hp[0], v = hp[1];
	ld tu = norm(u.HI - u.LO) * .5;
	ld mu = norm(u.HI + u.LO) * .5;
	if (!inside(u, mu)) mu = norm(mu + PI);
	ld du = r * cosl(tu);
	ld ru = r * sinl(tu);
	ld hu = r - du;
	ld tv = norm(v.HI - v.LO) * .5;
	ld mv = norm(v.HI + v.LO) * .5;
	if (!inside(v, mv)) mv = norm(mv + PI);
	ld dv = r * cosl(tv);
	ld rv = r * sinl(tv);
	ld hv = r - dv;
	if (eq(tu * 2, PI) && eq(tv * 2, PI)) {
		Polygon P;
		for (const Pos& p : hp) {
			if (p.LO < p.HI) P.push_back(p);
			else {
				P.push_back(Pos(0, p.LO));
				P.push_back(Pos(p.HI, 2 * PI));
			}
		}
		std::sort(P.begin(), P.end());
		P.push_back(Pos(2 * PI, 2 * PI));
		ld hi = 0, tt = 0;
		for (const Pos& p : P) {
			if (hi < p.LO) tt += (p.LO - hi), hi = p.HI;
			else hi = std::max(hi, p.HI);
		}
		ld vol = Sphere(0, 0, 0, r).vol();
		return vol * (tt / (2 * PI));
	}
	bool f1 = outer_check(u, v), f2 = outer_check(v, u);
	if (f1 && f2) {
		ld volu = Sphere(0, 0, 0, r).vol(hu);
		ld volv = Sphere(0, 0, 0, r).vol(hv);
		return Sphere(0, 0, 0, r).vol() - volu - volv;
	}
	if (f1) return Sphere(0, 0, 0, r).vol(r + r - hv);
	if (f2) return Sphere(0, 0, 0, r).vol(r + r - hu);
	Pos us = c.p(u.LO), ue = c.p(u.HI);
	Pos vs = c.p(v.LO), ve = c.p(v.HI);
	Seg U = Seg(us, ue);
	Seg V = Seg(vs, ve);
	Pos m = intersection(us, ue, vs, ve);
	ld dm = m.mag();
	ld tm = norm(m.rad());
	if (du < TOL && dv < TOL) { dm *= -1; tm = norm(tm + PI); }
	ld ttu = std::min(norm(tm - mu), norm(mu - tm));
	ld ttv = std::min(norm(tm - mv), norm(mv - tm));
	ld a_ = the(dm, r);
	ld suf = 0;
	ld x;
	x = intersection(U, V);
	assert(-TOL < x && x < 1 + TOL);
	x = .5 - x;
	if (inside(v, u.HI)) x *= -1;
	ld ang_u = the(x, 0.5);
	suf += Sphere(0, 0, 0, r).surf(hu) * ((PI * 2 - ang_u) / (PI * 2));
	x = intersection(V, U);
	assert(-TOL < x && x < 1 + TOL);
	x = .5 - x;
	if (inside(u, v.HI)) x *= -1;
	ld ang_v = the(x, 0.5);
	suf += Sphere(0, 0, 0, r).surf(hv) * ((PI * 2 - ang_v) / (PI * 2));
	int su = 0;
	if ((inside(Pos(u.LO, mu), tm) && inside(v, u.HI)) ||
		(inside(Pos(mu, u.HI), tm) && inside(v, u.LO))) su = -1;
	else su = 1;
	int sv = 0;
	if ((inside(Pos(v.LO, mv), tm) && inside(u, v.HI)) ||
		(inside(Pos(mv, v.HI), tm) && inside(u, v.LO))) sv = -1;
	else sv = 1;
	x = intersection(U, V);
	if (du < 0 && dv > 0 && (
		(inside(v, u.LO) && x > .5) ||
		(inside(v, u.HI) && x < .5)
		)) suf += (Sphere(0, 0, 0, r).surf() - area(a_, tu, tu, r, ttu)) * su;
	else suf += area(a_, tu, tu, r, ttu) * su;
	x = intersection(V, U);
	if (dv < 0 && du > 0 && (
		(inside(u, v.LO) && x > .5) ||
		(inside(u, v.HI) && x < .5)
		)) suf += (Sphere(0, 0, 0, r).surf() - area(a_, tv, tv, r, ttv)) * sv;
	else suf += area(a_, tv, tv, r, ttv) * sv;
	suf = Sphere(0, 0, 0, r).surf() - suf;
	ld ratio = suf / Sphere(0, 0, 0, r).surf();
	ld total = Sphere(0, 0, 0, r).vol() * ratio;
	total += cone_vol(ru, ang_u, du);
	total += cone_vol(rv, ang_v, dv);
	if (total < 0 || zero(du) || zero(dv)) {
		bool f = 1;
		Pos u_ = eq(tu * 2, PI) ? u : v;
		Pos v_ = eq(tu * 2, PI) ? v : u;
		assert(u_ != v_);
		ld suf = Sphere(0, 0, 0, r).surf() * .5, x;
		mv = norm(v_.HI + v_.LO) * .5;
		if (!inside(v_, mv)) mv = norm(mv + PI);
		if (eq(u_.LO, mv) || eq(u_.HI, mv)) {
			tv = norm(v_.HI - v_.LO) * .5;
			dv = r * cosl(tv);
			rv = r * sinl(tv);
			hv = r - dv;
			return Sphere(0, 0, 0, r).vol(r + r - hv) * .5;
		}
		if (inside(u_, mv)) { f = 0; std::swap(v_.LO, v_.HI); }
		mv = norm(v_.HI + v_.LO) * .5;
		if (!inside(v_, mv)) mv = norm(mv + PI);
		Pos us = c.p(u_.LO), ue = c.p(u_.HI);
		Pos vs = c.p(v_.LO), ve = c.p(v_.HI);
		Seg U = Seg(us, ue);
		Seg V = Seg(vs, ve);
		assert(!inside(u_, mv));
		tv = norm(v_.HI - v_.LO) * .5;
		dv = r * cosl(tv);
		rv = r * sinl(tv);
		hv = r - dv;
		x = intersection(U, V);
		assert(-TOL < x && x < 1 + TOL);
		x = .5 - x;
		if (inside(v_, u_.HI)) x *= -1;
		assert(x != 0);
		ld a_ = the(x, .5);
		x = intersection(V, U);
		assert(-TOL < x && x < 1 + TOL);
		x = .5 - x;
		if (inside(u_, v_.HI)) x *= -1;
		ld ang_v = the(x, 0.5);
		suf += Sphere(0, 0, 0, r).surf(hv) * ((PI * 2 - ang_v) / (PI * 2));
		ld ttv = 0;
		if (inside(u_, v_.HI)) ttv = std::min(norm(u_.LO - mv), norm(mv - u_.LO));
		else ttv = std::min(norm(u_.HI - mv), norm(mv - u_.HI));
		ld tri = area(a_, tv, tv, r, ttv);
		suf += tri;
		suf = Sphere(0, 0, 0, r).surf() - suf;
		ld ratio = suf / Sphere(0, 0, 0, r).surf();
		ld total = Sphere(0, 0, 0, r).vol() * ratio;
		total += cone_vol(rv, ang_v, dv);
		return f ? total : (Sphere(0, 0, 0, r).vol() * .5 - total);
	}
	return total;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	for (int i = 0; i < 3; i++)
		std::cin >> S[i].x >> S[i].y >> S[i].z >> S[i].r, F[i] = 0;
	std::sort(S, S + 3);
	ld dab = mag(S[0], S[1]);
	ld dac = mag(S[0], S[2]);
	Pos ca = Pos(0, 0);
	Pos cb = Pos(dab, 0);
	ld t = rad(S[0], S[1], S[2]);
	Pos cc = Pos(dac, 0).rot(t);
	C[0] = Circle(ca, S[0].r);
	C[1] = Circle(cb, S[1].r);
	C[2] = Circle(cc, S[2].r);
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				if (k == i || k == j) continue;
				Polygon B = tangent_convex_hull(C[i], C[j]);
				if (inner_check(B, C[k])) {
					std::cout << two_union(S[i], S[j]) << "\n";
					return;
				}
			}
		}
	}
	Vld T;
	for (int i = 0, i0, i1; i < 3; i++) {
		i0 = i;
		i1 = (i + 1) % 3;
		ld w = mag(S[i0], S[i1]);
		ld h = std::abs(S[i0].r - S[i1].r);
		T.push_back(sqrt(sq(w) + sq(h)));
	}
	Vld H = { (ld)S[0].r, (ld)S[1].r, (ld)S[2].r };
	ld A = tri_pil_vol(T, H);
	for (int i = 0; i < 3; i++) {
		int i0 = i, i1 = (i + 1) % 3, i2 = (i + 2) % 3;
		Pos hp1 = tangent(C[i0], C[i1]);
		Pos hp2 = tangent(C[i0], C[i2]);
		Polygon HP = { hp1, hp2 };
		if (!intersect(C[i0], hp1, hp2)) {
			A = 0;
			A += volume(S[i].r, HP);
			Pos hp3 = tangent(C[i1], C[i0]);
			A += volume(S[i1].r, { hp3 });
			Pos hp4 = tangent(C[i2], C[i0]);
			A += volume(S[i2].r, { hp4 });
			ld w0, w1, h0, w, h, t;
			Pos m0, m1;
			w0 = (C[i0].p(hp1.HI) - C[i0].p(hp1.LO)).mag();
			w1 = (C[i1].p(hp3.HI) - C[i1].p(hp3.LO)).mag();
			m0 = (C[i0].p(hp1.HI) + C[i0].p(hp1.LO)) * .5;
			m1 = (C[i1].p(hp3.HI) + C[i1].p(hp3.LO)) * .5;
			h0 = (m0 + m1).mag();
			w = std::abs(w0 - w1);
			t = atan2(h0, w);
			h = std::max(w0, w1) * tan(t);
			A += cone_vol(std::max(w0, w1), h, h0);
			w0 = (C[i0].p(hp2.HI) - C[i0].p(hp2.LO)).mag();
			w1 = (C[i2].p(hp4.HI) - C[i2].p(hp4.LO)).mag();
			m0 = (C[i0].p(hp2.HI) + C[i0].p(hp2.LO)) * .5;
			m1 = (C[i2].p(hp4.HI) + C[i2].p(hp4.LO)) * .5;
			h0 = (m0 + m1).mag();
			w = std::abs(w0 - w1);
			t = atan2(h0, w);
			h = std::max(w0, w1) * tan(t);
			A += cone_vol(std::max(w0, w1), h, h0);
			return;
		}
		A += volume(S[i].r, HP);
	}
	Pos3D c0 = pos3d(C[0].c);
	Pos3D p0 = pos3d(C[0]);
	Pos3D c1 = pos3d(C[1].c);
	Pos3D p1 = pos3d(C[1]);
	Pos3D c2 = pos3d(C[2].c);
	Pos3D p2 = pos3d(C[2]);
	Planar tn = planar(p2, p1, p0);
	Planar p01 = planar(p0, p1, c1);
	Planar p12 = planar(p1, p2, c2);
	Planar p20 = planar(p2, p0, c0);
	ld t01 = PI * 2 - rad(tn, p01) * 2;
	ld t12 = PI * 2 - rad(tn, p12) * 2;
	ld t20 = PI * 2 - rad(tn, p20) * 2;

	return;
}
int main() { solve(); return 0; }//boj23590
