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
inline int sign(const ld& x) { return x < 0 ? -1 : !!x; }
inline ld norm(ld th) {
	while (th < 0) th += 2 * PI;
	while (sign(th - 2 * PI) >= 0) th -= 2 * PI;
	return th;
}
inline int fit(const int& x, const int& lo, const int& hi) { return std::max(lo, std::min(hi, x)); }

#define ANG
//#define CAR
//#define RAD

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
Pos qry[105][3];
Vpos P[40];
bool cmp(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ll tq = p / q;
	return !tq ? p.Euc() <= q.Euc() : tq > 0;
}
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
std::vector<Pos> graham_scan(std::vector<Pos>& C) {
	std::vector<Pos> H;
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
bool convex_valiator(const std::vector<Pos>& H) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = H[i], & p1 = H[(i + 1) % sz], & p2 = H[(i + 2) % sz];
		std::cout << "p0:: " << p0 << " p1:: " << p1 << " p2:: " << p2 << "\n";
		int tq = ccw(p0, p1, p2);
		std::cout << "dir:: " << (tq > 0 ? "OK" : "BAD") << "\n";
		if (tq <= 0) return 0;
	}
	return 1;
}
int inner_check_bi_search(const std::vector<Pos>& H, const Pos& p) {//convex
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
bool two_polygon_check(const Vpos& U, const Vpos& V) {
	int szu = U.size();
	int szv = V.size();
	for (int i = 0; i < szu; i++) {
		const Pos& u1 = U[i], & u2 = U[(i + 1) % szu];
		if (inner_check_bi_search(V, u1) >= 0) return 0;
		for (int j = 0; j < szv; j++) {
			const Pos& v1 = V[j], & v2 = V[(j + 1) % szv];
			if (on_seg_strong(u1, u2, v1)) return 0;
			if (on_seg_strong(u1, u2, v2)) return 0;
			if (on_seg_strong(v1, v2, u1)) return 0;
			if (on_seg_strong(v1, v2, u2)) return 0;
		}
	}
	for (const Pos& v : V) if (inner_check_bi_search(U, v) >= 0) return 0;
	return 1;
}
bool val(const Pos& p, const int& i = -1) {
	bool f0 = -1000 <= p.x && p.x <= 1000 && -1000 <= p.y && p.y <= 1000;
	if (i == -1 || !f0) return f0;
	if (inner_check_bi_search(P[0], p) >= 0) return 0;
	return 1;
}
Pos fit(const Pos& p, const int& x) {
	return Pos(fit(p.x, -x, x), fit(p.y, -x, x));
}
Vpos Polygon_generator(const int& n, const int& i, const int& B, const Pos& cen, const int& tol) {
	int coo = 900 / B;
	Vpos C, H;
	Pos s = Pos(rnd.next(1, coo), 0);
	int cnt = 0;
	while (C.size() < n * 10) {
		cnt = cnt + 1;
		ld t = rnd.next() * 2 * PI;
		Pos p = cen + s.rot(t) + Pos(rnd.next(-tol, tol), rnd.next(-tol, tol));
		if (val(p, i) && std::find(C.begin(), C.end(), p) == C.end())
			C.push_back(p);
		if (cnt % 10000 == 0) {
			std::cout << "Now P[" << i << "] gene...\n";
			//std::cout << "C.sz:: " << C.size() << "\n";
		}
		if (cnt > 100000) return {};
	}
	H = graham_scan(C);
	//std::cout << "Random Hull gene\n";
	//std::cout << "H.sz:: bfr " << H.size() << "\n";
	if (H.size() > n) {
		//std::cout << "Hull is too big\n";
		int sz = H.size();
		int diff = sz - n;
		Vint I;
		while (I.size() < diff) {
			int i(rnd.next(0, sz - 1));
			if (std::find(I.begin(), I.end(), i) == I.end())
				I.push_back(i);
		}
		Vbool F(sz, 0);
		for (const int& i : I) F[i] = 1;
		Vpos V;
		for (int i = 0; i < sz; i++) if (!F[i]) V.push_back(H[i]);
		H = V;
	}
	return H;
}
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 1);
	//std::cout << "generating start\n";
	Vpos C, C1, C2, PR;
	int cnt = 0;
	C1 = { Pos(-20, -1000), Pos(-40, -999) };
	int x = -60;
	for (int i = 0; i < 48; i++) {
		int y = C1.back().y;
		C1.push_back(Pos(x, y));
		while (cross(C1[i], C1[i + 1], C1[i + 2]) >= 0) C1[i + 2].y++;
		if (i > 9) C1[i + 2].y++;
		x -= 20;
	}

	C2 = { Pos(20, -1000), Pos(40, -999) };
	x = 60;

	for (int i = 0; i < 48; i++) {
		int y = C2.back().y;
		C2.push_back(Pos(x, y));
		while (cross(C2[i], C2[i + 1], C2[i + 2]) <= 0) C2[i + 2].y++;
		if (i > 9) C2[i + 2].y++;
		x += 20;
	}

	for (Pos& p : C1) C.push_back(p);
	for (Pos& p : C2) C.push_back(p);
	//for (Pos& p : C) std::cout << p << "\n";
	Vpos H = graham_scan(C);
	std::cout << H.size() << "\n";
	std::cout << "100\n";
	for (Pos& p : C) (p).println();

	//Pos p0 = Pos(-999, 999);
	//Pos p1 = Pos(-980, 998);
	//Pos p2 = Pos(999, 999);
	//Pos v0 = Pos(1, 0);
	//Pos v1 = Pos(20, -1);
	//Pos v2 = Pos(-10, -1);
	//std::cout << "33\n";
	//for (int i = 0; i < 33; i++) {
	//	P[i] = { p0, p1, p2 };
	//	p0 = p1 + v0;
	//	p1 = p0 + v1;
	//	p2 = p2 + v2;
	//	std::cout << "3\n";
	//	for (Pos& p : P[i]) p.println();
	//}

	Pos p0 = Pos(-97, -100);
	Pos p1 = Pos(-96, -101);
	Pos p2 = Pos(-95, -100);
	std::cout << "33\n";
	for (int i = 0; i < 33; i++) {
		P[i] = { p0, p1, p2 };
		p0.x += + 6;
		p1.x += + 6;
		p2.x += + 6;
		std::cout << "3\n";
		for (Pos& p : P[i]) (p).println();
	}

	p0 = Pos(-350, -600);
	p1 = Pos(-50, -600);
	p2 = Pos(250, -600);
	std::cout << "100\n";
	for (int i = 0; i < 100; i++) {
		(p0).print();
		std::cout << " ";
		(p1).print();
		std::cout << " ";
		(p2).print();
		std::cout << "\n";
		p0.x += 1;
		p1.x += 1;
		p2.x += 1;
	}
	return 0;
}