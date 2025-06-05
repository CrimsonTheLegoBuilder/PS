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
#include <map>
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

int N;
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
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
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
struct Pos3D {
	ld x, y, z;
	int r, i;
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
std::vector<Pos3D> pos;
Pos3D S2C(const ld& lon, const ld& lat) {//Spherical to Cartesian
	ld phi = lon * PI / 180;
	ld the = lat * PI / 180;
	return Pos3D(cos(phi) * cos(the), sin(phi) * cos(the), sin(the));
}
Pos3D point(const Pos3D Xaxis, const Pos3D Yaxis, const ld& th) { return Xaxis * cos(th) + Yaxis * sin(th); }
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) { return norm(atan2(Yaxis * p, Xaxis * p)); }
//ld randTOL() {
//	ld rand01 = rand() / (ld)RAND_MAX;
//	ld err = (rand01 - .5) * TOL;
//	return err;
//}
//Pos3D add_noise(const Pos3D& p) {//refer to BIGINTEGER
//	return p + Pos3D(randTOL(), randTOL(), randTOL());
//}
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
//std::vector<Pos3D> graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
ld graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
	//if (C.size() < 3) {
	//	std::sort(C.begin(), C.end());
	//	return C;
	// }
	if (C.size() < 3) return 0;
	std::vector<Pos3D> H;
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pos3D& p, const Pos3D& q) -> bool {
		ld ret = ccw(C[0], p, q, norm);
		if (zero(ret)) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i], norm) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	//return H;
	return area(H, norm);
}
bool inner_check(const Pos3D& d1, const Pos3D& d2, const Pos3D& t) {
	Pos3D nrm = cross(O3D, d1, d2);
	Pos3D p1 = d1, p2 = d2;
	if (ccw(O3D, p1, p2, nrm) < 0) std::swap(p1, p2);
	return ccw(O3D, p1, t, nrm) >= 0 && ccw(O3D, p2, t, nrm) <= 0;
}
int inner_check(std::vector<Pos3D>& H, const Pos3D& p) {//for convex hull
	int sz = H.size();
	if (sz <= 1) return -1;
	if (sz == 2) {
		if (on_seg_strong(H[0], H[1], p)) return 0;
		else return -1;
	}
	Pos3D torque0 = cross(H[0], H[1], p);
	for (int i = 1; i < sz; i++) {
		Pos3D cur = H[i], nxt = H[(i + 1) % sz];
		Pos3D torqueI = cross(cur, nxt, p);
		if (zero(torqueI.mag())) {
			if (on_seg_strong(cur, nxt, p)) return 0;
			else return -1;
		}
		if (torque0 * torqueI < 0) return -1;
	}
	return 1;
}

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
int above(const Plane& S, const Pos3D& p) { return sign(p * S.norm() + S.d); }
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
bool circle_intersection(const Pos3D& a, const Pos3D& b, const ld& th, std::vector<Pos3D>& inxs) {
	inxs.clear();
	Pos3D mid = (a + b) * .5;
	if (zero(mid.mag())) return 0;
	ld x = cos(th) / mid.mag();
	if (x < -1 || 1 < x) return 0;
	Pos3D w = mid.unit() * x;
	ld ratio = sqrtl(1 - x * x);
	Pos3D h = (mid.unit() / (b - a).unit()) * ratio;
	inxs.push_back(w + h);
	inxs.push_back(w - h);
	return 1;
}
//bool plane_circle_intersection(const Pos3D& a, const Pos3D& perp, const ld& th, std::vector<Pos3D>& inxs) {
//	inxs.clear();
//	Pos3D vec = a - (perp * (perp * a));
//	if (zero(vec.mag())) return 0;
//	ld x = cos(th) / vec.mag();
//	if (x < -1 || 1 < x) return 0;
//	Pos3D w = vec.unit() * x;
//	ld ratio = sqrtl(1 - x * x);
//	Pos3D h = (vec.unit() / perp) * ratio;
//	inxs.push_back(w + h);
//	inxs.push_back(w - h);
//	return 1;
//}
int intersection(const Plane& p1, const Plane& p2, Line3D& l) {
	Pos3D n1 = p1.norm();
	Pos3D n2 = p2.norm();
	Pos3D dir = n2 / n1;
	if (zero(dir.mag())) {
		ld f = n1 * n2;
		ld d1 = dist(p1, O3D);
		ld d2 = dist(p2, O3D);
		if (sign(f) > 0) return sign(d2 - d1) >= 0 ? 0 : -1;
		else return sign(d2 + d1) >= 0 ? 0 : -2;
	}
	dir = dir.unit();
	Pos3D q1 = intersection(p1, Line3D(n1, O3D));
	Pos3D v1 = n1 / dir;
	Pos3D p0 = intersection(p2, Line3D(v1, q1));
	l = Line3D(dir, p0);
	return 1;
}
Pos3D s2c(const ld& phi, const ld& psi) {//Spherical to Cartesian
	ld lat = phi * PI / 180, lon = psi * PI / 180;
	return Pos3D(cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat));
}
bool plane_circle_intersection(const Pos3D& perp, const Pos3D& a, Polyhedron& inxs, const ld& R = 0) {
	inxs.clear();
	Pos3D vec = a - (perp * (perp * a));
	if (zero(vec.mag())) return 0;
	ld th = (ld)a.r / R;
	ld x = cos(th) / vec.mag();
	if (x < -1 || 1 < x) return 0;
	Pos3D w = vec.unit() * x;
	ld ratio = sqrtl(1 - x * x);
	Pos3D h = (vec.unit() / perp) * ratio;
	inxs.push_back(w + h);
	if (!zero(ratio)) inxs.push_back(w - h);
	return 1;
}
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) {
	ld X = Xaxis * p, Y = Yaxis * p;
	ld th = atan2(Y, X);
	return th;
}
ld angle(const Pos3D& a, const Pos3D& b) {
	Pos3D perp = (a / b).unit();
	Pos3D X = a.unit();//X-axis
	Pos3D Y = (perp / a).unit();//Y-axis
	return angle(X, Y, b);
}
bool inner_check(const Pos3D& d1, const Pos3D& d2, const Pos3D& t, const Pos3D& nrm) {
	return ccw(O3D, d1, t, nrm) >= 0 && ccw(O3D, d2, t, nrm) <= 0;
}
bool connectable(const Polyhedron& P, const Pos3D& a, const Pos3D& b, const int& i, const int& j, const ld& R = 0) {
	if (a == b) return 1;
	Pos3D perp = (a / b).unit();
	Polyhedron inxs;
	for (int k = 0; k < N; k++) {
		if (k == i || k == j) continue;
		Pos3D axis = (P[k] / perp).unit();
		if (zero(axis.Euc())) {
			if ((ld)P[k].r / R < PI * .5) continue;
			else return 0;
		}
		bool f = plane_circle_intersection(perp, P[k], inxs);
		if (!f) continue;
		if (inxs.size() == 1) {
			if (inner_check(a, b, inxs[0], perp)) return 0;
			continue;
		}
		Pos3D hi = inxs[0], lo = inxs[1];
		Pos3D mid = (perp / axis).unit();
		if (ccw(O3D, mid, lo, perp) < 0) std::swap(hi, lo);
		if (inner_check(lo, hi, a, perp) || inner_check(lo, hi, b, perp)) return 0;
		if (inner_check(a, b, lo, perp) || inner_check(a, b, hi, perp)) return 0;
	}
	return 1;
}
Polyhedron circle_circle_intersections(const Pos3D& p, const Pos3D& q, const ld& R = 0) {
	ld ang = angle(p, q);
	ld t1 = (ld)p.r / R;
	ld t2 = (ld)q.r / R;
	assert(t1 <= PI * .5); assert(t2 <= PI * .5);
	if (t1 + t2 < ang || std::abs(t1 - t2) >= ang) return {};
	ld d1 = cos(t1);
	ld d2 = cos(t2);
	Plane p1 = plane(p, -d1);
	Plane p2 = plane(q, -d2);
	Line3D inx;
	if (intersection(p1, p2, inx) != 1) return {};
	ld w = inx.p0.mag();
	ld h = sqrt(1 - sq(w));
	Pos3D v = inx.dir;
	return { inx.p0 + v * h, inx.p0 - v * h };
}
ld spherical_pythagorean(const ld& a, const ld& b, const ld& A, const ld& B) {
	assert(sin(a) > 0); assert(cos(b) > 0);
	ld cosc = fit(cos(a) / cos(b), -1, 1);
	ld c = acos(cosc);
	ld sinC = sin(c) * sin(A) / sin(a);
	return asin(sinC);
}
bool inner_check(const Pos3D& p, const Pos3D& q, const ld& R = 0) {
	ld ang = angle(p, q);
	ld t1 = (ld)p.r / R;
	ld t2 = (ld)q.r / R;
	assert(t1 <= PI * .5); assert(t2 <= PI * .5);
	if (p.r > q.r && std::abs(t1 - t2) >= ang) return 1;
	return 0;
}
int D_[13] = { 0, 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30 };
std::map<std::string, int> M;
void init() {
	M["jan"] = D_[2];
	M["feb"] = D_[2];
	M["mar"] = D_[3];
	M["apr"] = D_[4];
	M["may"] = D_[5];
	M["jun"] = D_[6];
	M["jul"] = D_[7];
	M["aug"] = D_[8];
	M["sep"] = D_[9];
	M["oct"] = D_[10];
	M["nov"] = D_[11];
	M["dec"] = D_[12];
	return;
}
ld w = 23.439281;
ld z = (2 * PI / 365) / 24;
ld rad(const ll& d, const std::string& s, const int& h) {
	return (M[s] + d) / 365.0 + z * h;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	ld l;
	std::cin >> l >> N;
	for (int i = 1; i < 13; i++) D_[i] = D_[i - 1] + D_[i];
	Pos3D e = Pos3D(acos(1));
	while (N--) {
		int D, h; std::string S;
		std::cin >> D >> S >> h;
		ld t = rad(D, S, h);


	}
}
int main() { solve(); return 0; }//boj31392
//boj 27712 10239 22635 29691 31392

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