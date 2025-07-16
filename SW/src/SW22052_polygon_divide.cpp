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
typedef long double ld;
//typedef double ld;
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

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

int N, M, Q;
ll A[LEN], T;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
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
ll area(Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) { if (area(H) < 0) std::reverse(H.begin(), H.end()); }
bool inner_check(Pos p0, Pos p1, Pos p2, const Pos& t) {
	if (ccw(p0, p1, p2) < 0) std::swap(p1, p2);
	return ccw(p0, p1, t) >= 0 && ccw(p1, p2, t) >= 0 && ccw(p2, p0, t) >= 0;
}
struct Pdd {
	ld x, y;
	Pdd(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	Pdd operator + (const Pdd& p) const { return { x + p.x, y + p.y }; }
	Pdd operator - (const Pdd& p) const { return { x - p.x, y - p.y }; }
	Pdd operator * (const ld& n) const { return { x * n, y * n }; }
	Pdd operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pdd& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pdd& p) const { return x * p.y - y * p.x; }
	Pdd operator ^ (const Pdd& p) const { return { x * p.x, y * p.y }; }
	Pdd& operator += (const Pdd& p) { x += p.x; y += p.y; return *this; }
	Pdd& operator -= (const Pdd& p) { x -= p.x; y -= p.y; return *this; }
	Pdd& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pdd& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	Pdd operator - () const { return { -x, -y }; }
	Pdd operator ~ () const { return { -y, x }; }
	Pdd operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pdd rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pdd unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pdd& p1, const Pdd& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pdd& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pdd& p) { os << p.x << " " << p.y; return os; }
};
typedef std::vector<Pdd> Vpdd;
Pdd conv(const Pos& p) { return Pdd(p.x, p.y); }
ld cross(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
Pdd intersection(const Pdd& p1, const Pdd& p2, const Pdd& q1, const Pdd& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld rad(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return rad(d2 - d1, d3 - d1); }
Pdd centroid(const Vpdd& H) {
	Pdd cen = Pdd(0, 0);
	ld a = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		ld a0 = H[i] / H[(i + 1) % sz];
		cen += (H[i] + H[(i + 1) % sz]) * a0;
		a += a0;
	}
	a *= .5;
	cen /= 6;
	if (!zero(a)) cen /= a;
	return cen;
}
ld ternary_search(const ll& a, const Pdd& I0, const Pdd& I1, const Pdd& J0, const Pdd& J1) {
	ld l = INF;
	int cnt = 50;
	Pdd v = I0 - I1;
	ld diff = A[N] * .5 - a;
	auto dist = [&](const ld& t) -> ld {
		Pdd i0 = I1 + v * t;
		ld tq0 = cross(i0, I1, J0);
		ld a0 = diff - tq0;
		ld d = (J0 - i0).mag();
		if (zero(d)) return INF;
		ld h = a0 / d;
		Pdd v1 = J0 - i0;
		Pdd i1 = i0 + ~v1.unit() * h;
		Pdd inx = intersection(J0, J1, i1, i1 + v1);
		if (on_seg_strong(J0, J1, inx)) return (i0 - inx).mag();
		return INF;
		};
	ld s = 0, e = 1;
	while (cnt--) {
		ld t1 = (s + s + e) / 3;
		ld t2 = (s + e + e) / 3;
		ld d1 = dist(t1);
		ld d2 = dist(t2);
		if (d2 > d1) e = t2;
		else s = t1;
	}
	return dist(s);
}
void query() {
	std::cin >> N; Polygon H(N); for (Pos& p : H) std::cin >> p;
	Vpdd P; for (Pos& p : H) P.push_back(conv(p));
	A[0] = 0;
	for (int i = 0; i < N; i++) A[i + 1] = A[i] + H[i] / H[(i + 1) % N];
	Pdd c = centroid(P);
	ld sh = INF, lg = -1;
	auto area_ = [&](const int& i, const int& j) -> ll {
		if (i <= j) return A[j] - A[i] + H[j] / H[i];
		ll a = A[i] - A[j] + H[i] / H[j];
		return A[N] - a;
		};
	for (int i = 0, j = 1; i < N; i++) {
		while (i != (j + 1) % N && (area_(i, (j + 1) % N) << 1) <= A[N])
			j = (j + 1) % N;
		int j1 = (j + 1) % N;

		ll a0 = area_(i, j);
		ll a1 = area_(j, i);
		assert(a0 <= a1);
		assert(a0 * 2 <= A[N]);

		ld diff = (A[N] * .5 - a0);
		Pdd I0 = conv(H[i]), J0 = conv(H[j]), J1 = conv(H[j1]);
		if (a0 == a1) {
			ld l0 = (I0 - J0).mag();
			sh = std::min(sh, l0);
			lg = std::max(lg, l0);
		}
		else {
			ll a = cross(H[i], H[j], H[j1]);
			if (a >= diff) {
				Pdd v = (J1 - J0);
				Pdd inx = J0 + (v * (diff / a));
				ld l0 = (inx - I0).mag();
				sh = std::min(sh, l0);
				lg = std::max(lg, l0);
			}
		}

		ll a = area_((i + 1) % N, j);

		Pdd I1 = conv(H[(i + 1) % N]);
		ld l1 = ternary_search(a, I0, I1, J0, J1);
		sh = std::min(sh, l1);

		I1 = conv(H[(i - 1 + N) % N]);
		if (A[N] != a0 * 2) {
			a = area_(i, j);
			ld l2 = ternary_search(a, I1, I0, J0, J1);
			sh = std::min(sh, l2);
		}
		else {
			J1 = J0;
			J0 = conv(H[(j - 1 + N) % N]);
			a = area_(i, (j - 1 + N) % N);
			ld l2 = ternary_search(a, I1, I0, J0, J1);
			sh = std::min(sh, l2);
		}
	}
	std::cout << sh << " " << lg << "\n";
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> Q; while (Q--) query();
	return;
}
int main() { solve(); return 0; }