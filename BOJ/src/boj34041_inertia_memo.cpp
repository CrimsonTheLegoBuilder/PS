#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <deque>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const int LEN = 1e5 + 10;
const ld INF = 1e17;
const ld TOL = 1e-7;
const ld PI = acos(-1);
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }

#define LINE 1
#define CIRCLE 2

#define LEFT 1
#define RIGHT -1

int N, Q;
ll A[LEN];
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
} C[LEN]; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld dist(const Pos& d1, const Pos& d2, const Pos& t, bool f = 0) {
	if (!f) return cross(d1, d2, t) / (d1 - d2).mag();
	if (sign(dot(d1, d2, t)) <= 0 &&
		sign(dot(d2, d1, t)) <= 0)
		return std::abs(cross(d1, d2, t)) / (d1 - d2).mag();
	return std::min((d1 - t).mag(), (d2 - t).mag());
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
	Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
struct Pii {
	int x, y;
	int i;
	Pii(int X = 0, int Y = 0) : x(X), y(Y) { i = -1; }
	bool operator == (const Pii& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pii& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pii& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pii operator + (const Pii& p) const { return { x + p.x, y + p.y }; }
	Pii operator - (const Pii& p) const { return { x - p.x, y - p.y }; }
	Pii operator * (const int& n) const { return { x * n, y * n }; }
	Pii operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pii& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pii& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	Pos p() const { return Pos(x, y); }
	friend std::istream& operator >> (std::istream& is, Pii& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pii& p) { os << p.x << " " << p.y; return os; }
};
typedef std::vector<Pii> Vpii;
ll cross(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return sign(cross(d1, d2, d3, d4)); }
bool collinear(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool on_seg_strong(const Pii& d1, const Pii& d2, const Pii& d3) { ll ret = dot(d1, d3, d2); return !ccw(d1, d2, d3) && ret >= 0; }
ld volume(const ld& a, const Pos& cen, const Seg& se) {
	ld d = std::abs(dist(se.s, se.e, cen));
	return 2 * PI * d * a;
}
ld query(const Vpii& P, const Polygon& H) {
	int qi, qj; std::cin >> qi >> qj; qi--; qj--;
	//std::cout << P[qi] << " " << P[qj] << "\n";
	//std::cout << collinear(P[qi], P[(qi + 1) % N], P[qj], P[(qj - 1 + N) % N]) << "\n";
	//std::cout << collinear(P[qj], P[(qj + 1) % N], P[qi], P[(qi - 1 + N) % N]) << "\n";
	if (collinear(P[qi], P[(qi + 1) % N], P[qj], P[(qj - 1 + N) % N])) return 0;
	if (collinear(P[qj], P[(qj + 1) % N], P[qi], P[(qi - 1 + N) % N])) return 0;
	if (qi == qj) return 0;
	if (qi > qj) std::swap(qi, qj);
	ll a = H[qj] / H[qi];
	ll A1 = A[qj] - A[qi] + a;
	ll A2 = A[N] - A1;
	if (!A1 || !A2) return 0;
	Pos c = (H[qi] + H[qj]);
	Pos C1 = C[qj] - C[qi] + c * a;
	Pos C2 = C[N] - (C[qj] - C[qi]) - c * a;
	C1 /= 6;
	C2 /= 6;
	ld a1 = A1 * .5;
	ld a2 = A2 * .5;
	C1 /= a1;
	C2 /= a2;
	Seg se = Seg(H[qi], H[qj]);
	return std::min(volume(a1, C1, se), volume(a2, C2, se));
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(20);
	std::cin >> N; Vpii P(N); for (Pii& p : P) std::cin >> p;
	Polygon H; for (Pii& p : P) H.push_back(p.p());
	A[0] = 0;
	C[0] = { 0, 0 };
	for (int i = 0; i < N; i++) {
		ll a = P[i] / P[(i + 1) % N];
		A[i + 1] = A[i] + a;
		C[i + 1] = C[i] + (H[i] + H[(i + 1) % N]) * a;
	}
	std::cin >> Q; while (Q--) {
		ld a = query(P, H);
		if (zero(a)) std::cout << "0\n";
		else std::cout << a << "\n";
	}
	return;
}
int main() { solve(); return 0; }