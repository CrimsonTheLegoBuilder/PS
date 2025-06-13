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
const int LEN = 1e5 + 1;
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
}; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
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
bool inner_check(const Polygon& H, const Pos& q) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) if (cross(H[i], H[(i + 1) % sz], q) <= 0) return 0;
	return 1;
}
bool inner_check(Pos p0, Pos p1, Pos p2, const Pos& q) {
	if (cross(p0, p1, p2) < 0) std::swap(p1, p2);
	return inner_check({ p0, p1, p2 }, q);
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
ll area(Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) { if (area(H) < 0) std::reverse(H.begin(), H.end()); }
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
};
struct Frac {
	ll num, den;//numerator, denominator (num/den)
	Frac(ll n = 0, ll d = 1) : num(n), den(d) {
		assert(den);
		if (den < 0) num *= -1, den *= -1;
		simplify();
	}
	void fit() {
		if (num < 0) num = 0, den = 1;
		else if (num > den) num = den = 1;
		return;
	}
	void simplify() {
		if (!num) { den = 1; return; }
		ll com = gcd(std::abs(num), std::abs(den));
		num /= com; den /= com;
		return;
	}
	friend std::ostream& operator << (std::ostream& os, const Frac& f) {
		os << f.num; if (f.den != 1) os << "/" << f.den; return os;
	}
	bool operator < (const Frac& o) const { return num * o.den < o.num * den; }
	bool operator <= (const Frac& o) const { return num * o.den <= o.num * den; }
	bool operator == (const Frac& o) const { return num == o.num && den == o.den; }
} _1 = Frac(1), _0 = Frac(0);
struct Range {
	Frac s, e;
	Range(Frac s_ = Frac(-1), Frac e_ = Frac(-1)) : s(s_), e(e_) {}
	bool operator < (const Range& r) const { return s == r.s ? e < r.e : s < r.s; }
};
typedef std::vector<Range> Vrange;
Frac intersection(const Seg& s1, const Seg& s2) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ll det = (q2 - q1) / (p2 - p1);
	if (!det) return Frac(-1);
	ll a1 = (q2 - q1) / (q1 - p1);
	Frac f1 = Frac(a1, det); f1.fit(); return f1;
	//ll a2 = (p2 - p1) / (p1 - q1);
	//Frac f2 = Frac(-a2, det);
	//if (_0 < f1 && f1 < _1 && _0 <= f2 && f2 <= _1) return f1;
	//return Frac(-1);
}
Range range(const Seg& s1, const Seg& s2, const Pos& o) {
	Pos p0 = o, p1 = s1.s, p2 = s1.e;
	Pos q1 = s2.s, q2 = s2.e;
	if (cross(p0, p1, p2) < 0) std::swap(p1, p2);
	if (cross(p0, q1, q2) < 0) std::swap(q1, q2);
	if (intersect(p0, p1, s2.s, s2.e) && intersect(p0, p2, s2.s, s2.e))
		return Range(_0, _1);
	bool i1 = inner_check(p0, p1, p2, q1);
	bool i2 = inner_check(p0, p1, p2, q2);
	if (!i1 && !i2) return Range();
	Frac f1 = _0, f2 = _1;
	if (i1) f1 = intersection(s1, Seg(p0, q1));
	if (i2) f2 = intersection(s1, Seg(p0, q2));
	return Range(f1, f2);
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	freopen("../../tests/aed/input/input0.txt", "r", stdin);
	freopen("../../tests/aed/input/ret.txt", "w", stdout);
	Pos u;
	std::cin >> N; Polygon H(N); std::cin >> u;
	for (Pos& p : H) std::cin >> p; norm(H);
	Vint I;
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
	std::cout << I.size() << "\n";
	for (int i : I) std::cout << i << " ";
	return;
}
int main() { solve(); return 0; }//boj29931
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