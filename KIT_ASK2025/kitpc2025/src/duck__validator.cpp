#define _CRT_SECURE_NO_WARNINGS
//#include "../inc/testlib.h"
#include "testlib.h"
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
const ld PI = acos(-1);
inline int sign(const ld& x) { return x < 0 ? -1 : !!x; }
inline ld norm(ld th) {
	while (th < 0) th += 2 * PI;
	while (sign(th - 2 * PI) >= 0) th -= 2 * PI;
	return th;
}
inline int fit(const int& x, const int& lo, const int& hi) { return std::max(lo, std::min(hi, x)); }
inline ll sq(const int& x) { return 1ll * x * x; }

int pt;
int R, L, N, M;
ll D;
struct Pos {
	int x, y, r, w;
	Pos(int x_ = 0, int y_ = 0, int r_ = 0, int w_ = 0) : x(x_), y(y_), r(r_), w(w_) {}
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
Vpos P, Q;
int C[4] = { 1000, 5000, 10000, 50000 };
bool meet(const Vpos& P, const int& i, const int& j, const int& l = 0) {
	ll d = (P[i] - P[j]).Euc();
	ll r1 = sq(P[i].r + P[j].r + l);
	return d <= r1;
}
int main(int argc, char* argv[]) {
	registerValidation(argc, argv);
	int x = 0, y = 0, r = 0, w = 0;

	R = inf.readInt(1, 1'000'000'000, "R");
	inf.readSpace();
	L = inf.readInt(1, 1'000'000'000, "L");
	inf.skipBlanks();
	N = inf.readInt(0, 100, "N");
	inf.readSpace();
	M = inf.readInt(0, 100, "M");
	inf.skipBlanks();

	for (int i = 0; i < N; i++) {
		x = inf.readInt(-1'000'000'000, 1'000'000'000, "x");
		inf.readSpace();
		y = inf.readInt(-1'000'000'000, 1'000'000'000, "y");
		inf.readSpace();
		r = inf.readInt(1, 1'000'000'000, "r");
		inf.skipBlanks();
		Pos p = Pos(x, y, r);
		bool f = p.Euc() <= sq(R);
		ensuref(f, "A few leaves are outside the pond\n");
		P.push_back(Pos(x, y, r));
	}

	for (int j = 1; j <= M; j++) {
		x = inf.readInt(-1'000'000'000, 1'000'000'000, "x");
		inf.readSpace();
		y = inf.readInt(-1'000'000'000, 1'000'000'000, "y");
		inf.readSpace();
		w = inf.readInt(1000, 50000, "r");
		inf.skipBlanks();
		bool f = 0;
		for (int k = 0; k < 4; k++) if (w == C[k]) f = 1;
		ensuref(f, "Currency Unit Error\n");
		Q.push_back(Pos(x, y, 0, w));
	}

	inf.readEof();
	return 0;
}