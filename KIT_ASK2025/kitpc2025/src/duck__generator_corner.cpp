#define _CRT_SECURE_NO_WARNINGS
#include "../inc/testlib.h"
//#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
typedef long long ll;
typedef long double ld;
typedef std::vector<bool> Vbool;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::vector<ld> Vld;
const int LEN = 105;
const ld PI = acos(-1);
inline int sign(const ld& x) { return x < 0 ? -1 : !!x; }
inline ld norm(ld th) {
	while (th < 0) th += 2 * PI;
	while (sign(th - 2 * PI) >= 0) th -= 2 * PI;
	return th;
}
inline int fit(const int& x, const int& lo, const int& hi) { return std::max(lo, std::min(hi, x)); }
inline ll sq(const ll& x) { return 1ll * x * x; }

int pt;
int R, L, N, M;
ll D;
struct Pos {
	int x, y, r, w;
	Pos(int x_ = 0, int y_ = 0, int r_ = 0, int w_ = 0) : x(x_), y(y_), r(r_) { w = w_; }
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
	Pos rot(const ld& the) const { return { int(x * cos(the) - y * sin(the)), int(x * sin(the) + y * cos(the)) }; }
	int quad() const { return y > 0 || y == 0 && x >= 0; }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
	void println() const { std::cout << x << " " << y << "\n";  return; }
	void print() const { std::cout << x << " " << y;  return; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Vpos;
Pos P[LEN], Q[LEN];
int C[4] = { 1000, 5000, 10000, 50000 };
#define RND 20
int main(int argc, char* argv[]) {
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(20);
	ld l = 0;
	l = sqrtl(sq(14144006) + sq(14140265));
	std::cout << l << "\n";
	l = sqrtl(sq(141429548) + sq(141413164));
	std::cout << l << "\n";
	l = sqrtl(sq(707123023) + sq(707090539));
	std::cout << l << "\n";
	l = sqrtl(sq(50000000) - sq(20));
	std::cout << l << "\n";
	l = sqrtl(sq(49999999) + sq(20));
	std::cout << l << "\n";
	l = sqrtl(sq(99999999) - sq(1));
	//9999999800000001 1
	//9,999,999,800,000,000
	//9,999,999,800,000,001
	std::cout << l << "\n";
	std::cout << sq(99999999) << " " << sq(1) << "\n";

	//registerGen(argc, argv, 0);
	std::cout << "generating start\n";
	for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
	R = 1e9;
	L = 1;
	N = 0;
	M = 0;
	D = sq(R);

	P[N++] = Pos(0, -999999999, 99999998);
	P[N++] = Pos(1, -900000000, 1);
	P[N++] = Pos(0, 999999999, 99999998);
	P[N++] = Pos(1, 900000000, 1);
	P[N++] = Pos(999999999, 0, 99999998);
	P[N++] = Pos(900000000, 1, 1);
	P[N++] = Pos(-999999999, 0, 99999998);
	P[N++] = Pos(-900000000, 1, 1);

	Q[M++] = Pos(1, -899999999, 0, 10000);
	Q[M++] = Pos(1, 899999999, 0, 10000);
	Q[M++] = Pos(899999999, 1, 0, 10000);
	Q[M++] = Pos(-899999999, 1, 0, 10000);
	Q[M++] = Pos(707106780, 707106780, 0, 10000);
	Q[M++] = Pos(-707106780, 707106780, 0, 10000);
	Q[M++] = Pos(-707106780, -707106780, 0, 10000);
	Q[M++] = Pos(707106780, -707106780, 0, 10000);

	std::cout << "generating done\n";
	std::cout << "\n\n\n";
	freopen("../tests/duck/in/25.in", "w", stdout);
	std::cout << R << " " << L << "\n";
	std::cout << N << " " << M << "\n";
	for (int i = 0; i < N; i++) std::cout << P[i] << " " << P[i].r << "\n";
	for (int i = 0; i < M; i++) std::cout << Q[i] << " " << Q[i].w << "\n";
	return 0;
}