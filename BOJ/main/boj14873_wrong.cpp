#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <deque>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
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

int N, M, K, T, Q;
struct Pos {
	ld x, y;
	Pos(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& n) const { return { x * n, y * n }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos rot(const ld& t) const { return Pos(x * cosl(t) - y * sinl(t), x * sinl(t) + y * cosl(t)); }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrtl(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2l(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
}; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
ld dist(const Pos& d1, const Pos& d2, const Pos& t, bool f = 0) {
	if (!f) return cross(d1, d2, t) / (d1 - d2).mag();
	if (sign(projection(d1, d2, d2, t)) <= 0 &&
		sign(projection(d2, d1, d1, t)) <= 0)
		return std::abs(cross(d1, d2, t)) / (d1 - d2).mag();
	return std::min((d1 - t).mag(), (d2 - t).mag());
}
ld area(const Polygon& H) {
	ld A = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) A += H[i] / H[(i + 1) % sz];
	return A * .5;
}
void norm(Polygon& H) {
	ld A = area(H); if (A < 0) std::reverse(H.begin(), H.end());
	return;
}
struct Seg {
	Pos s, e, dir;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) { dir = e - s; }
	Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
typedef std::vector<Seg> Vseg;
struct Circle {
	Pos c;
	ld r;
	Circle(Pos c_ = Pos(), ld r_ = 0) : c(c_), r(r_) {}
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld rad(const Pos& p) const { return (p - c).rad(); }
	ld area(const ld& lo = 0, const ld& hi = 2 * PI) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		//if (hi < lo) { return green(lo, 2 * PI) + green(0, hi); }
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
	bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
	bool operator < (const Pos& p) const { return r < (c - p).mag(); }
};
typedef std::vector<Circle> Disks;
Vld intersections(const Circle& a, const Circle& b) {
	Pos ca = a.c, cb = b.c;
	Pos vec = cb - ca;
	ld ra = a.r, rb = b.r;
	if (zero(a.r) || zero(b.r)) return {};
	ld distance = vec.mag();
	ld rd = vec.rad();
	if (vec.Euc() > sq(ra + rb) + TOL) return {};
	if (vec.Euc() < sq(ra - rb) - TOL) return {};
	ld X = (ra * ra - rb * rb + vec.Euc()) / (2 * distance * ra);
	if (X < -1) X = -1;
	if (X > 1) X = 1;
	ld h = acos(X);
	Vld ret = {};
	ret.push_back(norm(rd - h));
	if (zero(h)) return ret;
	ret.push_back(norm(rd + h));
	return ret;
}
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
		if (0 < hi && hi < 1) ret.push_back(hi);
		if (zero(det)) return ret;
		if (0 < lo && lo < 1) ret.push_back(lo);
	}
	else {
		auto the = [&](ld rt) { return q.rad(s + (e - s) * rt); };
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
		if (zero(det)) return ret;
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
	}
	return ret;
}
ld green(const Circle& c, const ld& lo, const ld& hi) {
	if (lo < hi) return c.green(lo, hi);
	return c.green(lo, 2 * PI) + c.green(0, hi);
}
struct Arc {
	Circle c;
	ld lo, hi;
	Arc(Circle c_, ld l_, ld h_) : c(c_), lo(l_), hi(h_) {}
	bool inside(const Pos& p) const {
		ld t = c.rad(p);
		if (lo < hi) return lo < t && t < hi;
		else return lo < t || t < hi;
	}
};
typedef std::vector<Arc> Arcs;
Arcs AR, AL;
ld area(const Arc& q) { return q.c.area(0, norm(q.hi - q.lo)); }
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(5);
	Circle C; std::cin >> C.c.x >> C.c.y >> C.r;
	std::cin >> N; Polygon H(N); for (Pos& p : H) std::cin >> p.x >> p.y; norm(H);
	int r = 0, l = 0;
	bool meet = 0;
	bool f0 = 0;
	ld A = 0;
	for (int i = 0, j; i < N; i++) {
		if (ccw(C.c, H[r], H[i]) < 0 || (!ccw(C.c, H[r], H[i]) && dot(C.c, H[r], H[i]) < 0)) r = i;
		if (ccw(C.c, H[l], H[i]) > 0 || (!ccw(C.c, H[l], H[i]) && dot(C.c, H[l], H[i]) < 0)) l = i;
		if (dist(H[i], H[(i + 1) % N], C.c, 1) < C.r) meet = 1;
		j = (i + 1) % N;
		Vld inxs = circle_line_intersections(C, H[j], H[i], CIRCLE);
		if (inxs.size() == 2) {
			ld lo = norm(inxs[0]), hi = norm(inxs[1]);
			ld t = norm(hi - lo);
			A += C.area(0, t);
			A += std::abs(C.r * C.r * .5 * sin(t));
			std::cout << A << "\n";
			return;
		}
	}
	if (!meet) { std::cout << C.area() << "\n"; return; }
	ld lo = norm((H[l] - C.c).rad());
	for (int j = l, k; 1; j = (j + 1) % N) {
		k = (j + 1) % N;
		bool fj = C >= H[j], fk = C >= H[k];
		if (fj) break;
		if (fk) {
			Vld inx = circle_line_intersections(C, H[k], H[j], CIRCLE);
			assert(inx.size() == 1);
			lo = inx[0];
			l = j;
			break;
		}
	}
	ld hi = norm((H[r] - C.c).rad());
	for (int j = r, i; 1; j = (j - 1 + N) % N) {
		i = (j - 1 + N) % N;
		bool fj = C >= H[j], fi = C >= H[i];
		if (fj) break;
		if (fi) {
			Vld inx = circle_line_intersections(C, H[i], H[j], CIRCLE);
			assert(inx.size() == 1);
			hi = inx[0];
			r = j;
			break;
		}
	}
	A = 0;
	for (int i = l, j; i != r; i = (i + 1) % N) {
		j = (i + 1) % N;
		bool fi = C >= H[i], fj = C >= H[j];
		if (fi && fj) {
			A += cross(C.c, H[j], H[i]) * .5;
			continue;
		}
		Vld inx = circle_line_intersections(C, H[j], H[i]);
		assert(inx.size() == 1);
		ld x = inx[0];
		if (fi) A += cross(C.c, Seg(H[j], H[i]).p(x), H[i]) * .5;
		if (fj) A += cross(C.c, H[j], Seg(H[j], H[i]).p(x)) * .5;
	}
	A += C.area(0, norm(hi - lo));
	ld L = C.r;
	int el = l, er = r;
	Circle cl, cr;
	ld aa = 0;
	for (int j = l, i, k; 1; j = (j - 1 + N) % N) {
		i = (j - 1 + N) % N;
		k = (j + 1) % N;
		lo = norm((H[i] - H[j]).rad());
		hi = norm((H[j] - H[k]).rad());
		if (j == l) hi = norm((H[j] - C.c).rad());
		ld d = (H[j] - H[k]).mag();
		if (j == l) d = (H[j] - C.c).mag();
		if (L <= d) break;
		L -= d;
		Circle c = Circle(H[j], L);
		AL.push_back(Arc(c, lo, hi));
		if (c < H[i]) { el = j; break; }
	}
	ld R = C.r;
	for (int j = r, i, k; 1; j = (j + 1) % N) {
		i = (j - 1 + N) % N;
		k = (j + 1) % N;
		lo = norm((H[j] - H[i]).rad());
		hi = norm((H[k] - H[j]).rad());
		if (j == r) lo = norm((H[j] - C.c).rad());
		ld d = (H[j] - H[i]).mag();
		if (j == r) d = (H[j] - C.c).mag();
		if (R <= d) break;
		R -= d;
		Circle c = Circle(H[j], R);
		AR.push_back(Arc(c, lo, hi));
		if (c < H[k]) { er = j; break; }
	}
	//for (Circle& c : ) {
	//	std::cout << c.c.x << c.c.y << c.r << "\n";
	//}
	//std::cout << "FUCK::\n";
	Vld inxs = intersections(AR.back().c, AL.back().c);
	ld x = 0;
	//std::cout << "FUCK::\n";
	if ((er + 1) % N == el && inxs.size() == 2) {
		std::deque<Pos> dq;
		while (1) {
			int szl = AL.size();
			int szr = AR.size();
			if (szl < 2 || szr < 2) {
				std::cout << "something wrong::\n";
				return;
			}
			inxs = intersections(AR.back().c, AL.back().c);
			x = inxs[0];
			Pos m = AR.back().c.p(x);
			if (!AR.back().inside(m)) { AR.pop_back(); dq.push_front(m); continue; }
			if (!AL.back().inside(m)) { AL.pop_back(); dq.push_back(m); continue; }
		}
	}
	for (const Arc& q : AL) A += area(q);
	for (const Arc& q : AR) A += area(q);
	std::cout << A << "\n";
	return;
}
int main() { solve(); return 0; }//boj14873
//boj30123 27712 3607 10239

/*

20 90 328
10
40 90
80 210
130 230
190 220
210 180
210 100
180 50
130 30
80 20
60 40

150 160 100
4
10 60
150 120
300 50
150 10

1 6 18
5
6 -1
3 0
2 6
3 9
6 13

10 60 180
5
60 -10
30 0
20 60
30 90
60 130
90858.58632

7 15 32
11
13 16
12 20
12 24
14 27
17 29
22 31
27 31
29 27
27 19
21 14
17 13

0 0 4
3
0 2
1 4
-1 4

*/
