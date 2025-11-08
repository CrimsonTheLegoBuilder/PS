#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <deque>
#include <map>
#include <set>
#include <random>
#include <array>
#include <tuple>
#include <complex>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::pair<int, int> pi;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-7;
const ld PI = acos(-1);
const int LEN = 105;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define __FUCK__ ;
#define WHAT_THE_FUCK
//#define DEBUG

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2
#define SEG 3
#define POS 4

#define STRONG 0
#define WEAK 1

int N, P, R;
struct Pos {
	ld x, y;
	int i, f;
	bool rv;
	Pos(ld x_ = 0, ld y_ = 0, int f_ = -1) : x(x_), y(y_), f(f_) { i = -1, rv = 0; }
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	//bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	//bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	//bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return *this < p || *this == p; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& scalar) const { return { x * scalar, y * scalar }; }
	Pos operator / (const ld& scalar) const { return { x / scalar, y / scalar }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& scale) { x *= scale; y *= scale; return *this; }
	Pos& operator /= (const ld& scale) { x /= scale; y /= scale; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pos rot(ld the) { return { x * cos(the) - y * sin(the), x * sin(the) + y * cos(the) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2l(y, x); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "(" << p.x << ", " << p.y << ")"; return os; }
} p0, p1, key, vec; const Pos O = Pos(0, 0); const Pos INVAL = Pos(INF, INF);
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = cross(d1, d2, d3); return sign(ret); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = dot(d1, d3, d2); return !ccw(d1, d2, d3) && sign(ret) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = dot(d1, d3, d2); return !ccw(d1, d2, d3) && sign(ret) > 0; }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
struct Seg {
	Pos a, b;
	int i, f;
	Seg(Pos a_ = Pos(), Pos b_ = Pos(), int i_ = -1, int f_ = -1) : a(a_), b(b_), i(i_), f(f_) {}
	Pos inx(const Seg& o) const { return intersection(a, b, o.a, o.b); }
	Pos p(const ld& rt = .5) const { return a + (b - a) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (a.x - b.x);
	}
};
struct Circle {
	Pos c;
	ld r;
	Circle(Pos c_ = Pos(), ld r_ = 0) : c(c_), r(r_) {}
	bool operator == (const Circle& q) const { return c == q.c && r == q.r; }
	bool operator != (const Circle& q) const { return !(*this == q); }
	bool operator < (const Circle& q) const { return c == q.c ? r < q.r : c < q.c; }
	//bool operator < (const Circle& q) const { return r < q.r && (c - q.c).mag() + r < q.r + TOL; }
	bool outside(const Circle& q) const { return sign((c - q.c).Euc() - sq(r + q.r)) >= 0; }
	Circle operator + (const Circle& q) const { return { c + q.c, r + q.r }; }
	Circle operator - (const Circle& q) const { return { c - q.c, r - q.r }; }
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld rad(const Pos& p) const { return norm((p - c).rad()); }
	ld area(const ld& lo = 0, const ld& hi = 2 * PI) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
	ld H(const ld& th) const { return sin(th) * c.x + cos(th) * c.y + r; }//coord trans | check right
	//bool operator < (const Pos& p) const { return r < (c - p).mag(); }
	bool operator < (const Pos& p) const { return sign(r - (c - p).mag()) < 0; }
	bool operator > (const Pos& p) const { return r > (c - p).mag(); }
	bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
	friend std::istream& operator >> (std::istream& is, Circle& c) { is >> c.c >> c.r; return is; }
	friend std::ostream& operator << (std::ostream& os, const Circle& c) { os << c.c << " " << c.r; return os; }
};
Vld circle_line_intersections(const Circle& q, const Pos& s, const Pos& e, const int& f = LINE) {
	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
	Pos vec = e - s;
	Pos OM = s - q.c;
	ld a = vec.Euc();
	ld b = vec * OM;
	ld c = OM.Euc() - q.r * q.r;
	ld J = b * b - a * c;
	if (J < -TOL) return {};
	ld det = sqrt(std::max((ld)0, J));
	ld lo = (-b - det) / a;
	ld hi = (-b + det) / a;
	Vld ret;
	if (f == LINE) {
		ret.push_back(hi);
		//if (-TOL < hi && hi < 1 + TOL) ret.push_back(hi);
		if (zero(det)) return ret;
		ret.push_back(lo);
		//if (-TOL < lo && lo < 1 + TOL) ret.push_back(lo);
	}
	else {
		auto the = [&](ld rt) { return norm(q.rad(s + (e - s) * rt)); };
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
		if (zero(det)) return ret;
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
	}
	return ret;
}
Pos shadow(const Circle& pp, const Circle& rd, const Circle& og) {
	Pos v = og.c - rd.c;
	ld a = v.mag(), b = og.r + rd.r;
	ld t = acos(b / a);

	Vld inxs;
	ld hi, lo;

	Pos r1 = v.rot(t) / a * rd.r + rd.c;
	Pos p1 = (-v).rot(t) / a * og.r + og.c;
	inxs = circle_line_intersections(pp, r1, p1, LINE);
	assert(inxs.size() == 2);
	hi = inxs[0]; lo = inxs[1];
	if (hi < lo) std::swap(hi, lo);
	assert(hi > 0 && lo < 0);
	Pos x1 = Seg(r1, p1).p(hi);

	Pos r2 = v.rot(-t) / a * rd.r + rd.c;
	Pos p2 = (-v).rot(-t) / a * og.r + og.c;
	inxs = circle_line_intersections(pp, r2, p2, LINE);
	assert(inxs.size() == 2);
	hi = inxs[0]; lo = inxs[1];
	if (hi < lo) std::swap(hi, lo);
	assert(hi > 0 && lo < 0);
	Pos x2 = Seg(r2, p2).p(hi);

	ld t1 = pp.rad(x1);
	ld t2 = pp.rad(x2);
	return Pos(t1, t2);
}
Circle C[LEN];
bool query() {
	int x, y, r;
	std::cin >> N >> P >> x >> y >> r;
	if (!N && !P && !x && !y && !r) return 0;
	Circle pp = Circle(Pos(0, 0), P), rd = Circle(Pos(x, y), r);
	for (int i = 0; i < N; i++) std::cin >> C[i];
	Polygon vp = { Pos(0, 0) };
	for (int i = 0; i < N; i++) {
		Pos sd = shadow(pp, rd, C[i]);
		//std::cout << "shadow:: " << sd << "\n";
		ld hi = sd.HI, lo = sd.LO;
		if (lo < hi) vp.push_back(sd);
		else {
			vp.push_back(Pos(lo, PI * 2));
			vp.push_back(Pos(0, hi));
		}
	}
	vp.push_back(Pos(2 * PI, 2 * PI));
	std::sort(vp.begin(), vp.end());
	ld hi = 0;
	ld t = 0;
	for (const Pos& p : vp) {
		if (hi < p.LO) t += (p.LO - hi), hi = p.HI;
		else hi = std::max(hi, p.HI);
	}
	std::cout << t / (PI * 2) << "\n";
	return 1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(4);
	while (query());
	return;
}
int main() { solve(); return 0; }//boj4793
