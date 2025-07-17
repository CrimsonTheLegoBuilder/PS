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
const ld TOL = 1e-15;
const ld PI = acos(-1);
const int LEN = 1 << 7;
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
#define SEG 3
#define POS 4

#define STRONG 0
#define WEAK 1

int N, M, K, T, Q;
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
Polygon P[LEN];
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
struct Seg {
	Pos s, e, dir;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) { dir = e - s; }
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
ld intersection(const Seg& s1, const Seg& s2, const bool& f = 0) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == 1) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
Pos get_pos(const Pos& l, const Seg& p, const Seg& q) {
	Pos p1 = p.s, p2 = p.e, q1 = q.s, q2 = q.e;
	if (intersect(l, q1, p1, p2)) return Pos(0, 0);
	if (intersect(l, q2, p1, p2)) return Pos(0, 0);
	if (!ccw(p1, p2, l)) {
		if (sign(on_seg_strong(p1, l, p2))) std::swap(p1, p2);
		if (intersect(p1, l, q1, q2, WEAK)) return Pos(0, 1);
		return Pos(0, 0);
	}
	if (ccw(p2, l, p1) < 0) std::swap(p1, p2);
	if (!inside(p2, l, p1, q1, WEAK) && !inside(p2, l, p1, q2, WEAK)) {
		if (intersect(l, p1, q1, q2) && intersect(l, p2, q1, q2)) return Pos(0, 1);
		else return Pos(0, 0);
	}
	Polygon tri = { p1, p2, l };
	ld r1 = 0, r2 = 0;
	bool in1 = inside(q2, l, q1, p.s), in2 = inside(q2, l, q1, p.e);
	if (!in1 && !in2) {
		r1 = intersection(p, Seg(l, q1), 1);
		r2 = intersection(p, Seg(l, q2), 1);
		if (r2 < r1) std::swap(r1, r2);
	}
	else if (in1 && in2) return r1 = 0, r2 = 1;
	else if (in1) {
		r1 = 0;
		if (inner_check(tri, q1)) r2 = intersection(p, Seg(l, q1), 1);
		if (inner_check(tri, q2)) r2 = intersection(p, Seg(l, q2), 1);
	}
	else if (in2) {
		r2 = 1;
		if (inner_check(tri, q1)) r1 = intersection(p, Seg(l, q1), 1);
		if (inner_check(tri, q2)) r1 = intersection(p, Seg(l, q2), 1);
	}
	if (r2 < r1) std::swap(r1, r2);
	return Pos(r1, r2);
}
bool query() {
	std::cin >> N;
	if (!N) return 0;
	ld t = 0;
	Pos s, e; std::cin >> s >> e;
	Seg se = Seg(s, e);
	for (int n = 0; n < N; n++) {
		P[n].resize(4);
		for (Pos& p : P[n]) std::cin >> p;
		norm(P[n]);
	}
	const Polygon& P0 = P[0];
	Polygon R;
	for (int n = 1; n < N; n++) {
		int ir = -1, il = -1;
		int jr = -1, jl = -1;
		const Polygon& H = P[n];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				int i0 = (i + 3) % 4, i1 = i, i2 = (i + 1) % 4;
				int j0 = (j + 3) % 4, j1 = j, j2 = (j + 1) % 4;
				if (ccw(H[j1], P0[i1], P0[i0]) <= 0
					&& ccw(H[j1], P0[i1], P0[i2]) <= 0
					&& ccw(P0[i1], H[j1], H[j0]) <= 0
					&& ccw(P0[i1], H[j1], H[j2]) <= 0) {
					il = i1;
					jr = j1;
				}
				if (ccw(H[j1], P0[i1], P0[i0]) >= 0
					&& ccw(H[j1], P0[i1], P0[i2]) >= 0
					&& ccw(P0[i1], H[j1], H[j0]) >= 0
					&& ccw(P0[i1], H[j1], H[j2]) >= 0) {
					ir = i1;
					jl = j1;
				}
			}
		}
		assert(ir != -1);
		assert(il != -1);
		assert(jr != -1);
		assert(jl != -1);
		Pos m = intersection(P0[ir], H[jl], H[jr], P0[il]);
		Pos JR = H[jr];
		Pos JL = H[jl];
		if (jr == jl) {
			JR = m + (m - P0[il]).unit();
			JL = m + (m - P0[ir]).unit();
		}
		Pos inx = get_pos(m, se, Seg(JL, JR));
		if (!eq(inx.x, inx.y)) R.push_back(inx);
	}
	ld hi = 0;
	ld ratio = 0;
	std::sort(R.begin(), R.end());
	R.push_back(Pos(1, 1));
	for (const Pos& r : R) {
		if (r.LO > hi) ratio += (r.LO - hi), hi = r.HI;
		else hi = std::max(hi, r.HI);
	}
	std::cout << (s - e).mag() * ratio << "\n";
	return 1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	freopen("../../tests/train.in", "r", stdin);
	freopen("../../tests/train.txt", "w", stdout);
	while (query());
	return;
}
int main() { solve(); return 0; }//11465

/*
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
const ld TOL = 1e-15;
const ld PI = acos(-1);
const int LEN = 1 << 7;
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

int N, M, K, T, Q;
struct Pos {
	ld x, y;
	int h, i;
	Pos(ld x_ = 0, ld y_ = 0, int h_ = 0, int i_ = 0) : x(x_), y(y_), h(h_), i(i_) {}
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
Polygon P[LEN];
Pos cen[LEN];
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool cmpt(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ld tq = p / q;
	return !tq ? p.Euc() < q.Euc() : tq > 0;
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
ld intersection(const Seg& s1, const Seg& s2, const bool& f = 0) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == 1) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
Pos get_pos(const Pos& l, const Seg& p, const Seg& q) {
	Pos p1 = p.s, p2 = p.e, q1 = q.s, q2 = q.e;
	if (intersect(l, q1, p1, p2)) return Pos(0, 0);
	if (intersect(l, q2, p1, p2)) return Pos(0, 0);
	if (!ccw(p1, p2, l)) {
		if (sign(on_seg_strong(p1, l, p2))) std::swap(p1, p2);
		if (intersect(p1, l, q1, q2, WEAK)) return Pos(0, 1);
		return Pos(0, 0);
	}
	if (ccw(p2, l, p1) < 0) std::swap(p1, p2);
	//if (ccw(l, p1, p2) < 0) std::swap(p1, p2);
	if (!inside(p2, l, p1, q1, WEAK) && !inside(p2, l, p1, q2, WEAK)) {
		//std::cout << "FUCK::\n";
		//std::cout << "l:: " << l << "\n";
		//std::cout << "p1:: " << p1 << "\n";
		//std::cout << "p2:: " << p2 << "\n";
		//std::cout << "q1:: " << q1 << "\n";
		//std::cout << "q2:: " << q2 << "\n";
		if (intersect(l, p1, q1, q2) && intersect(l, p2, q1, q2)) return Pos(0, 1);
		else return Pos(0, 0);
	}
	Polygon tri = { p1, p2, l };
	ld r1 = 0, r2 = 0;
	bool in1 = inside(q2, l, q1, p.s), in2 = inside(q2, l, q1, p.e);
	if (!in1 && !in2) {
		r1 = intersection(p, Seg(l, q1), 1);
		r2 = intersection(p, Seg(l, q2), 1);
		if (r2 < r1) std::swap(r1, r2);
	}
	else if (in1 && in2) return r1 = 0, r2 = 1;
	else if (in1) {
		r1 = 0;
		if (inner_check(tri, q1)) r2 = intersection(p, Seg(l, q1), 1);
		if (inner_check(tri, q2)) r2 = intersection(p, Seg(l, q2), 1);
	}
	else if (in2) {
		r2 = 1;
		if (inner_check(tri, q1)) r1 = intersection(p, Seg(l, q1), 1);
		if (inner_check(tri, q2)) r1 = intersection(p, Seg(l, q2), 1);
	}
	if (r2 < r1) std::swap(r1, r2);
	return Pos(r1, r2);
}
bool query() {
	//Q++;
	//std::cout << Q << "\n";
	std::cin >> N;
	if (!N) return 0;
	ld t = 0;
	Pos s, e; std::cin >> s >> e;
	//std::cout << s << " " << e << "\n";
	Seg se = Seg(s, e);
	for (int n = 0; n < N; n++) {
		P[n].resize(4);
		for (Pos& p : P[n]) std::cin >> p;
		//for (Pos& p : P[n]) std::cout << p << " ";
		//std::cout << "\n";
		norm(P[n]);
		for (int i = 0; i < 4; i++) P[n][i].h = n, P[n][i].i = i;
		cen[n] = (P[n][0] + P[n][2]) * .5;
	}
	const Polygon& P0 = P[0];
	Polygon R;
	for (int n = 1; n < N; n++) {
		int ir = -1, il = -1;
		int jr = -1, jl = -1;
		const Polygon& H = P[n];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				int i0 = (i + 3) % 4, i1 = i, i2 = (i + 1) % 4;
				int j0 = (j + 3) % 4, j1 = j, j2 = (j + 1) % 4;
				if (ccw(H[j1], P0[i1], P0[i0]) <= 0
					&& ccw(H[j1], P0[i1], P0[i2]) <= 0
					&& ccw(P0[i1], H[j1], H[j0]) <= 0
					&& ccw(P0[i1], H[j1], H[j2]) <= 0) {
					il = i1;
					jr = j1;
				}
				if (ccw(H[j1], P0[i1], P0[i0]) >= 0
					&& ccw(H[j1], P0[i1], P0[i2]) >= 0
					&& ccw(P0[i1], H[j1], H[j0]) >= 0
					&& ccw(P0[i1], H[j1], H[j2]) >= 0) {
					ir = i1;
					jl = j1;
				}
			}
		}
		assert(ir != -1);
		assert(il != -1);
		assert(jr != -1);
		assert(jl != -1);
		//std::cout << "ir:: " << ir << "\n";
		//std::cout << "il:: " << il << "\n";
		//std::cout << "jr:: " << jr << "\n";
		//std::cout << "jl:: " << jl << "\n";
		//std::cout << "P0[ir]:: " << P0[ir] << "\n";
		//std::cout << "P0[il]:: " << P0[il] << "\n";
		//std::cout << "H[jr]:: " << H[jr] << "\n";
		//std::cout << "H[jl]:: " << H[jl] << "\n";
		Pos m = intersection(P0[ir], H[jl], H[jr], P0[il]);
		Pos JR = H[jr];
		Pos JL = H[jl];
		if (jr == jl) {
			JR = m + (m - P0[il]).unit();
			JL = m + (m - P0[ir]).unit();
		}
		Pos inx = get_pos(m, se, Seg(JL, JR));
		//std::cout << "inx[" << n << "]:: " << inx << "\n";
		if (!eq(inx.x, inx.y)) R.push_back(inx);
	}
	ld hi = 0;
	ld ratio = 0;
	//for (Pos& r : R) std::cout << "r:: " << r << "\n";
	std::sort(R.begin(), R.end());
	R.push_back(Pos(1, 1));
	for (const Pos& r : R) {
		if (r.LO > hi) ratio += (r.LO - hi), hi = r.HI;
		else hi = std::max(hi, r.HI); 
	}
	std::cout << (s - e).mag() * ratio << "\n";
	return 1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	freopen("../../tests/train.in", "r", stdin);
	freopen("../../tests/train.txt", "w", stdout);
	while (query());
	return;
}
int main() { solve(); return 0; }//11465
//boj 27712 10239 22635 29691 31392 16068
*/

/*

2
0 0 0 1
3 -4 11 -8 13 -4 5 0
6 -6 6 -9 -6 -9 -6 -6
0

2
10000. 10000. 10000. -10000.
-10000. 10000. -10000. 9999. -9999. 9999. -9999. 10000.
9999. -9999. 9999. -10000. 9998. -10000. 9998. -9999.
0

2
10000. -10000. -10000. -10000.
10000. 10000. 9999. 10000. 9999. 9999. 10000. 9999.
-9999. -9999. -10000. -9999. -10000. -9998. -9999. -9998.
0

20000.

1
8
2 2
-2 -2
2 -2
-2 2
3 0
-3 0
0 3
0 -3

1
10
-2 -4
5 2
7 7
10 -2
9 -5
-5 10
1 -5
4 -9
5 1
7 -7

1
4
4 5
-2 4
1 4
-5 -2

1
5
1 1
2 4
1 4
0 4
-1 1

3
4
0 3
1 3
3 1
3 0
4
-4 0
5 3
0 -4
-1 0
5
4 4
5 0
3 3
3 2
-4 2
5
1 1
2 4
1 4
0 4
-1 1

2
4
-1 -1
1 1
1 -1
-1 1
5
-1 -1
1 1
1 -1
-1 1
0 2


*/

#define BOJ
#ifdef BOJ
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
//int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068