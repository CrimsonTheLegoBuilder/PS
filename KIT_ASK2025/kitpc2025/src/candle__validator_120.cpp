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
const ld PI = acos(-1);
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline ld norm(ld th) {
	while (th < 0) th += 2 * PI;
	while (sign(th - 2 * PI) >= 0) th -= 2 * PI;
	return th;
}
inline int fit(const int& x, const int& lo, const int& hi) { return std::max(lo, std::min(hi, x)); }

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
	Pos rot(const ld& the) const { return { int(x * cos(the) - y * sin(the)), int(x * sin(the) + y * cos(the)) }; }
	int quad() const { return y > 0 || y == 0 && x >= 0; }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "(" << p.x << ", " << p.y << ")"; return os; }
	void println() const { std::cout << x << " " << y << "\n";  return; }
	void print() const { std::cout << x << " " << y;  return; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Vpos;
Pos qry[200][3];
Vpos P[50];
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
bool convex_valiator(const std::vector<Pos>& H) {
	int sz = H.size();
	ll a = 0;
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = H[i], & p1 = H[(i + 1) % sz], & p2 = H[(i + 2) % sz];
		ll tq = cross(p0, p1, p2);
		if (tq <= 0) return 0;
	}
	return 1;
}
int inner_check(const std::vector<Pos>& H, const Pos& p) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		int i0 = i, i1 = (i + 1) % sz;
		if (on_seg_strong(H[i0], H[i1], p)) return 1;
		if (ccw(H[i0], H[i1], p) < 0) return 0;
	}
	return 2;
}
bool two_polygon_check(const Vpos& U, const Vpos& V) {
	int szu = U.size();
	int szv = V.size();
	for (int i = 0; i < szu; i++) {
		const Pos& u1 = U[i], & u2 = U[(i + 1) % szu];
		if (inner_check(V, u1) > 0) return 0;
		for (int j = 0; j < szv; j++) {
			const Pos& v1 = V[j], & v2 = V[(j + 1) % szv];
			if (on_seg_strong(u1, u2, v1)) return 0;
			if (on_seg_strong(u1, u2, v2)) return 0;
			if (on_seg_strong(v1, v2, u1)) return 0;
			if (on_seg_strong(v1, v2, u2)) return 0;
			if (intersect(u1, u2, v1, v2)) return 0;
		}
	}
	for (const Pos& v : V) if (inner_check(U, v) > 0) return 0;
	return 1;
}
int main(int argc, char* argv[]) {
	registerValidation(argc, argv);

	int x = 0, y = 0;

	int n = inf.readInt(3, 120, "N");
	inf.skipBlanks();

	for (int i = 0; i < n; i++) {
		x = inf.readInt(-1000, 1000, "x");
		inf.readSpace();
		y = inf.readInt(-1000, 1000, "y");
		inf.skipBlanks();
		P[0].push_back(Pos(x, y));
	}

	int k = inf.readInt(1, 40, "K");
	inf.skipBlanks();

	int sm = 0;
	for (int j = 1; j <= k; j++) {
		int m = inf.readInt(3, 120, "M");
		inf.skipBlanks();
		sm += m;
		for (int l = 0; l < m; l++) {
			x = inf.readInt(-1000, 1000, "x");
			inf.readSpace();
			y = inf.readInt(-1000, 1000, "y");
			inf.skipBlanks();
			P[j].push_back(Pos(x, y));
		}
	}
	ensuref(sm <= 120, "Summation of all m is too big!! ( > 100)", &sm);

	int q = inf.readInt(1, 120, "Q");
	inf.skipBlanks();

	for (int i = 0; i < q; i++) {
		x = inf.readInt(-1000, 1000, "rx");
		inf.readSpace();
		y = inf.readInt(-1000, 1000, "ry");
		inf.readSpace();
		qry[i][0] = Pos(x, y);

		x = inf.readInt(-1000, 1000, "gx");
		inf.readSpace();
		y = inf.readInt(-1000, 1000, "gy");
		inf.readSpace();
		qry[i][1] = Pos(x, y);

		x = inf.readInt(-1000, 1000, "bx");
		inf.readSpace();
		y = inf.readInt(-1000, 1000, "by");
		inf.skipBlanks();
		qry[i][2] = Pos(x, y);
	}

	for (int i = 0; i <= k; i++) {
		bool f = convex_valiator(P[i]);
		ensuref(f, "Polygon %d is not convex!!", i);
	}

	for (int i = 1; i <= k; i++) {
		for (Pos& p : P[i]) {
			int f = inner_check(P[0], p);
			ensuref(f == 2, "Blue polygon %d is not strictly inside of RED CONVEX POLYGON!!", i);
		}
		for (int j = i + 1; j <= k; j++) {
			bool f = two_polygon_check(P[i], P[j]);
			ensuref(f, "Blue polygon %d and %d are overlap!!", i, j);
		}
	}

	for (int i = 0; i < q; i++) {
		for (int j = 0; j < 3; j++) {
			Pos p = qry[i][j];
			for (int l = 0; l <= k; l++) {
				int sz = P[l].size();
				for (int o = 0; o < sz; o++) {
					int o0 = o, o1 = (o + 1) % sz;
					bool f = on_seg_strong(P[l][o0], P[l][o1], p);
					ensuref(!f, "candle %d - %d is on segment!!", i, j);
				}
			}
		}
	}

	inf.readEof();
	return 0;
}