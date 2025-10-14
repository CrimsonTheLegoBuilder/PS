#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <random>
#include <array>
#include <tuple>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 3e4 + 5;
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
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
//ll pow_fuck(ll a, ll b) {
//	ll ret = 1;
//	while (b) {
//		if (b & 1) ret = ret * a % MOD;
//		a = a * a % MOD;
//		b >>= 1;
//	}
//	return ret;
//}
//ll powmod(ll a, ll b) {
//	ll res = 1; a %= MOD;
//	assert(b >= 0);
//	for (; b; b >>= 1) {
//		if (b & 1) res = res * a % MOD;
//		a = a * a % MOD;
//	}
//	return res;
//}
struct Info { ll area, l, r; };

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

//freopen("../../../input_data/triathlon_tests/triath.20", "r", stdin);
//freopen("../../../input_data/triathlon_tests/triathlon_out.txt", "w", stdout);

//Euler characteristic : v - e + f == 1
//Pick`s Theorem : A = i + b/2 - 1

//2D============================================================================//

int N, M, T, Q, R;
struct Pos {
	int x, y;
	//ll x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
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
} S[LEN], E[LEN]; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool cmpt(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ll tq = p / q;
	return !tq ? p.Euc() < q.Euc() : tq > 0;
}
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
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
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
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
ld intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) {
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	return ((q2 - q1) / (q1 - p1)) / det;
}
bool inner_check(const Polygon& H, const Pos& q, const Pos& d = Pos(0, 0)) {//convex
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos& p1 = H[i], & p2 = H[(i + 1) % sz];
		if (ccw(p1, p2, q) < 0) return 0;
		if (on_seg_strong(p1, p2, q) && !eq(d.x, d.y)) {
			if (sign((p2 - p1) / d) > 0) return 1;
			else return 0;
		}
		if (on_seg_strong(p1, p2, q)) return 2;
	}
	return 1;
}
Pos shadow(const Pos& l, const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2, const int& n) {
	if (!inside(p2, l, p1, q1, WEAK) && !inside(p2, l, p1, q2, WEAK)) {
		if (intersect(l, p1, q1, q2) && intersect(l, p2, q1, q2)) return Pos(0, n);
		else return Pos(-1, -1);
	}
	Polygon tri = { p1, p2, l };
	bool in1 = inner_check(tri, q1), in2 = inner_check(tri, q2);
	if (!in1 && !in2) return Pos(-1, -1);
	ld r1 = 0, r2 = n;
	if (in1 && in2) {
		r1 = intersection(p1, p1, l, q1);
		r1 *= n;
		if ()
		r2 = intersection(p1, p2, l, q2);
	}
	else if (in1) {
		r1 = intersection(p1, p2, l, q1);
	}
	else if (in2) {
		r2 = intersection(p1, p2, l, q2);
	}
	else r1 = r2 = -1;
	if (r2 < r1) std::swap(r1, r2);
	return Pos(r1, r2);
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	Pos J;
	std::cin >> N >> J >> R;
	for (int i = 0; i < R; i++) {
		int n; std::cin >> n;
		Polygon P(n);
		for (Pos& p : P) std::cin >> p;
		Pos s = P[0], e = P[0];
		for (const Pos& p : P) {
			if (ccw(J, s, p) < 0) s = p;
			if (ccw(J, e, p) > 0) e = p;
		}
		S[i] = s; E[i] = e;
	}
	Polygon B = { Pos(0, 0), Pos(N, 0), Pos(N, N), Pos(0, N) };
	for (int b = 0; b < 4; b++) {
		Pos s = B[b], e = B[(b + 1) % 4];
		for (int i = 0; i < R; i++) {
			Pos a = 
		}
	}
}
int main() { solve(); return 0; }