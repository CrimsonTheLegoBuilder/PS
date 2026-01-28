#define _CRT_SECURE_NO_WARNINGS 
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-6;
const ld PI = acos(-1);
const int LEN = 10005;
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
#define STRONG 0
#define WEAK 1
#define LO x
#define HI y
#define LINE 1
#define CIRCLE 2

ld heron(const ld& a, const ld& b, const ld& c) {
	ld s = (a + b + c) / 2;
	ld v = s * (s - a) * (s - b) * (s - c);
	v = std::max(v, (ld)0);
	ld ret = sqrt(v);
	return ret;
}
ld rad(const ld& r, const ld& d) {
	if (zero(r - d)) return PI;
	ld dif = r - d;
	bool f = dif < 0;
	if (f) dif *= -1;
	ld h = sqrt(r * r - dif * dif);
	ld t = atan2(h, dif);
	if (f) t = PI - t;
	return t * 2;
}
ld vol(const ll& r) { return (4. / 3) * PI * r * r * r; }
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
	Pos rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
bool inner_check(const Polygon& H, const Pos& p) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) if (ccw(H[i], H[(i + 1) % sz], p) < 0) return 0;
	return 1;
}
ld area(const Polygon& H) {
	ld a = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a * .5;
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
};
Vld circle_line_intersections(const Circle& q, const Pos& s, const Pos& e, const int& f = LINE) {
	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
	Pos vec = e - s;
	Pos OM = s - q.c;
	ld a = vec.Euc();
	ld b = vec * OM;
	ld c = OM.Euc() - q.r * q.r;
	ld J = b * b - a * c;
	//std::cout << "J:::: " << J << "\n";
	if (zero(J)) return { .5 };
	if (J < -TOL) return {};
	ld det = sqrt(std::max((ld)0, J));
	ld lo = (-b - det) / a;
	ld hi = (-b + det) / a;
	Vld ret;
	if (f == LINE) {
		if (0 < hi && hi < 1) ret.push_back(hi);
		if (zero(det)) return ret;
		if (0 < lo && lo < 1) ret.push_back(lo);
	}
	else {
		auto the = [&](ld rt) { return norm(q.rad(s + (e - s) * rt)); };
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
		if (zero(det)) return ret;
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
	}
	return ret;
}
struct Sphere {
	ll x, y, z, r;
	Sphere(ll x_ = 0, ll y_ = 0, ll z_ = 0, ll r_ = 0) : x(x_), y(y_), z(z_), r(r_) {}
	bool operator < (const Sphere& q) const { return x == q.x ? y == q.y ? z < q.z : y < q.y : x < q.x; }
	Sphere operator - (const Sphere& q) const { return { x - q.x, y - q.y, z - q.z, 0 }; }
	ll operator * (const Sphere& q) const { return x * q.x + y * q.y + z * q.z; }
	ld vol() const { return (4. / 3) * PI * r * r * r; }
	ld vol(const ld& h) const { return PI * h * h * (3 * r - h) / 3; }
	ld surf() const { return PI * 4 * r * r; }
	ld surf(const ld& h) const { return PI * 2 * r * h; }
} S[3];
bool collinear() {
	ll x1 = S[1].x - S[0].x;
	ll y1 = S[1].y - S[0].y;
	ll z1 = S[1].z - S[0].z;
	ll x2 = S[2].x - S[1].x;
	ll y2 = S[2].y - S[1].y;
	ll z2 = S[2].z - S[1].z;
	return (y1 * z2 - z1 * y2) == 0 && (z1 * x2 - x1 * z2) == 0 && (x1 * y2 - y1 * x2) == 0;
}
ll Euc(const Sphere& p, const Sphere& q) { return sq(p.x - q.x) + sq(p.y - q.y) + sq(p.z - q.z); }
ld mag(const Sphere& p, const Sphere& q) { return sqrtl(Euc(p, q)); }
ld rad(const Sphere& a, const Sphere& b, const Sphere& c) {
	ld dab = mag(a, b);
	ld dac = mag(a, c);
	ld det = (b - a) * (c - a);
	ld proj = det / dab;
	ld ret = fit(proj / dac, -1, 1);
	return acos(ret);
}
Pos centroid(const Polygon& H) {
	Pos cen = Pos(0, 0);
	ld a = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		ld tq = H[i] / H[(i + 1) % sz];
		cen += (H[i] + H[(i + 1) % sz]) * tq;
		a += tq;
	}
	a *= .5;
	cen /= 6;
	if (!zero(a)) cen /= a;
	return cen;
}
ld rot_vol(const Polygon& H, const Pos& s1, const Pos& s2, const ld& t = PI * 2) {
	Pos cen = centroid(H);
	ld d = std::abs(cross(s1, s2, cen) / (s1 - s2).mag());
	ld a = area(H);
	return a * d * t;
}
bool inner_check(const Circle& c, Pos ca, Pos p1, Pos cb, Pos q1) {
	if (ccw(ca, p1, c.c) >= 0 || ccw(cb, q1, c.c) <= 0) return 0;
	if (ccw(q1, p1, c.c) <= 0) return 0;
	ld d = cross(q1, p1, c.c) / (q1 - p1).mag();
	return sign(d - c.r) >= 0;
}
void get_tangent(const Circle& c1, const Circle& c2, Pos& p1, Pos& p2, Pos& q1, Pos& q2) {
	Pos vec = c2.c - c1.c;
	ll w = c1.r - c2.r;
	if (!w) {
		Pos v = ~vec.unit() * c1.r;
		p1 = c1.c + v;
		p2 = c1.c - v;
		q1 = c2.c + v;
		q2 = c2.c - v;
		return;
	}
	bool f = w < 0;
	if (f) w *= -1;
	ld d = vec.mag();
	ld h = sqrt(d * d - w * w);
	ld t = atan2(h, w);
	if (f) t = PI - t;
	Pos v;
	v = vec.unit() * c1.r;
	p1 = c1.c + v.rot(t);
	p2 = c1.c + v.rot(-t);
	v = vec.unit() * c2.r;
	q1 = c2.c + v.rot(t);
	q2 = c2.c + v.rot(-t);
	return;
}
ld two_convex_hull(const Circle& c1, const Circle& c2) {
	//std::cout << "FUCK::\n";
	if (eq(c1.r, c2.r)) return vol(c1.r) + c1.r * c1.r * PI * (c1.c - c2.c).mag();
	Pos p1, p2, q1, q2;
	get_tangent(c1, c2, p1, p2, q1, q2);
	Pos m1 = (p1 + p2) * .5;
	Pos m2 = (q1 + q2) * .5;
	Polygon H = { p1, m1, m2, q1 };
	ld V = rot_vol(H, m1, m2) + vol(c1.r) + vol(c2.r);
	ld h, d;
	d = cross(p2, p1, c1.c) / (p1 - p2).mag();
	h = c1.r - d;
	V -= Sphere(0, 0, 0, c1.r).vol(h);
	d = cross(q1, q2, c2.c) / (q2 - q1).mag();
	h = c2.r - d;
	V -= Sphere(0, 0, 0, c2.r).vol(h);
	return V;
}
void spherical_triangle_angles(const ld& a, const ld& b, const ld& c, ld& A_, ld& B_, ld& C_) {
	A_ = acos(fit((cos(a) - cos(b) * cos(c)) / (sin(b) * sin(c)), -1, 1));
	B_ = acos(fit((cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c)), -1, 1));
	C_ = acos(fit((cos(c) - cos(a) * cos(b)) / (sin(a) * sin(b)), -1, 1));
	return;
}
//ld area(const ld& a, const ld& b, const ld& c, const ll& r, const ld& t = 0) {
ld area(const ld& a, const ld& b, const ld& c, const ll& r) {
	ld A_, B_, C_;
	//if (a >= PI) {
	//	spherical_triangle_angles(a * .5, t, c, A_, B_, C_);
	//	return r * r * (A_ + B_ + C_ - PI) * 2;
	//}
	spherical_triangle_angles(a, b, c, A_, B_, C_);
	return r * r * (A_ + B_ + C_ - PI);
}
ld corner_vol(const ll& r, const ld& h1, const ld& h2, const ld& t1, const ld& t2, const ld& r1, const ld& r2, const ld& r3) {
	Sphere s = Sphere(0, 0, 0, r);
	ld suf = s.surf();
	ld s1 = s.surf(h1) * (t1 / (2 * PI));
	ld s2 = s.surf(h2) * (t2 / (2 * PI));
	ld s3 = area(r1, r2, r3, r) * 2;
	//std::cout << "DEBUG::\n r:: " << r << " suf:: " << suf << " s1:: " << s1 << " s2:: " << s2 << " s3 " << s3 << " s.surf():: " << s.surf() << "\n";
	suf -= s1;
	suf -= s2;
	suf -= s3;
	ld rat = suf / s.surf();
	return s.vol() * rat;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	for (int i = 0; i < 3; i++) std::cin >> S[i].x >> S[i].y >> S[i].z >> S[i].r;
	std::sort(S, S + 3);
	if (collinear()) {
		ld dab = mag(S[0], S[2]);
		ld dac = mag(S[0], S[1]);
		Pos ca = Pos(0, 0);
		Pos cb = Pos(dab, 0);
		Pos cc = Pos(dac, 0);
		Circle c1 = Circle(ca, S[0].r);
		Circle c2 = Circle(cb, S[2].r);
		Circle c3 = Circle(cc, S[1].r);
		Pos p1, p2, q1, q2;
		get_tangent(c1, c2, p1, p2, q1, q2);
		ld V = 0;
		if (inner_check(c3, ca, p1, cb, q1)) V = two_convex_hull(c1, c2);
		else V = two_convex_hull(c1, c3) + two_convex_hull(c2, c3) - vol(c3.r);
		std::cout << V << "\n";
		return;
	}
	for (int i = 0, j, k; i < 3; i++) {
		j = (i + 1) % 3;
		k = (j + 1) % 3;
		ld dab = mag(S[i], S[j]);
		ld dac = mag(S[i], S[k]);
		Pos ca = Pos(0, 0);
		Pos cb = Pos(dab, 0);
		ld t = rad(S[i], S[j], S[k]);
		Pos cc = Pos(dac, 0).rot(t);
		Circle c1 = Circle(ca, S[i].r);
		Circle c2 = Circle(cb, S[j].r);
		Circle c3 = Circle(cc, S[k].r);
		Pos p1, p2, q1, q2;
		get_tangent(c1, c2, p1, p2, q1, q2);
		if (inner_check(c3, ca, p1, cb, q1)) {
			std::cout << two_convex_hull(c1, c2) << "\n";
			return;
		}
		//std::cout << "c1:: " << c1 << "\n";
		//std::cout << "c2:: " << c2 << "\n";
		//std::cout << "c3:: " << c3 << "\n";
		//std::cout << "p1:: " << p1  << " p2:: " << p2 << "\n";
		//std::cout << "q1:: " << q1  << " q2:: " << q2 << "\n";
		Vld t1 = circle_line_intersections(c3, p1, q1);
		Vld t2 = circle_line_intersections(c3, p2, q2);
		//std::cout << "t1.size():: " << t1.size() << "\n";
		//std::cout << "t2.size():: " << t2.size() << "\n";
		if (t1.size() && t2.size()) {;
			ld V = two_convex_hull(c1, c3) + two_convex_hull(c2, c3) - vol(c3.r);
			std::cout << V << "\n";
			return;
		}
	}
	ld V = 0;
	ld h = (S[0].r + S[1].r + S[2].r) / 3.;
	ld d1 = sqrt(Euc(S[0], S[1]) - sq(S[0].r - S[1].r));
	ld d2 = sqrt(Euc(S[1], S[2]) - sq(S[1].r - S[2].r));
	ld d3 = sqrt(Euc(S[2], S[0]) - sq(S[2].r - S[0].r));
	V += heron(d1, d2, d3) * h * 2;
	//std::cout << "DEBUG:: " << h << "\n";
	//std::cout << "DEBUG:: " << V << "\n";
	for (int i = 0, j, k; i < 3; i++) {
		j = (i + 1) % 3;
		k = (j + 1) % 3;
		ld dab = mag(S[i], S[j]);
		ld dac = mag(S[i], S[k]);
		Pos ca = Pos(0, 0);
		Pos cb = Pos(dab, 0);
		ld t = rad(S[i], S[j], S[k]);
		Pos cc = Pos(dac, 0).rot(t);
		Circle c1 = Circle(ca, S[i].r);
		Circle c2 = Circle(cb, S[j].r);
		Circle c3 = Circle(cc, S[k].r);
		Pos p1, p2, q1, q2;
		get_tangent(c1, c2, p1, p2, q1, q2);
		Pos u1, u2, w1, w2;
		get_tangent(c1, c3, u1, u2, w1, w2);
		Pos inx = intersection(p1, p2, u1, u2);

		ld r, l, h1, h2, t1, t2, r1, r2, r3;
		r3 = t;
		
		l = (p1 - p2).mag();
		r = l * .5;
		h = cross(p2, p1, ca) / l;
		h1 = S[i].r - h;
		t1 = rad(r, (inx - p2).mag());

		Pos a1, a2, b1, b2;
		get_tangent(c2, c3, a1, a2, b1, b2);
		Pos xx = intersection(a1, a2, q1, q2);
		l = (q1 - q2).mag();
		r = l * .5;
		//std::cout << "DEBUG:: t1:: " << t1 << " t11:: " << rad(r, (xx - q2).mag()) <<"\n";
		r1 = rad(cb - ca, p1 - ca);

		Polygon H = { ca, p2, q2, cb };
		V += rot_vol(H, ca, cb, t1);
		//std::cout << "DEBUG:: " << V << "\n";

		l = (u1 - u2).mag();
		r = l * .5;
		h = cross(u2, u1, ca) / l;
		h2 = S[i].r - h;
		t2 = rad(r, (inx - u1).mag());
		r2 = rad(cc - ca, u1 - ca);

		V += corner_vol(c1.r, h1, h2, t1, t2, r1, r2, r3);
		//std::cout << "DEBUG:: " << V << "\n";
	}
	std::cout << V << "\n";
	return;
}
int main() { solve(); return 0; }//boj23590

/*

0 1 0 1
30 10 0 10
60 2 0 2

-100 -99 -100 1
0 0 0 100
100 -98 100 2

*/
