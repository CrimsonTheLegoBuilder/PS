#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-9;
const ld PI = acosl(-1);
const int LEN = 505;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo = 0, const ld& hi = 1) { return std::min(hi, std::max(lo, x)); }

#define STRONG 0
#define WEAK 1

#define LINE 1
#define CIRCLE 2

#define LO x
#define HI y

int N, V;
ld A[LEN][LEN];
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
bool inner_check(const Polygon& H, const Pos& p) {
	int sz = H.size();
	for (int i = 0; i < sz; i++)
		if (ccw(H[i], H[(i + 1) % sz], p) <= 0)
			return 0;
	return 1;
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
struct Circle {
	Pos c;
	ld r;
	Circle(Pos c_ = Pos(), ld r_ = 0) : c(c_), r(r_) {}
	bool operator == (const Circle& C) const { return c == C.c && eq(r, C.r); }
	//bool operator >= (const Pos& p) const { return sign(r * r - (c - p).Euc()) >= 0; }
	bool operator >= (const Pos& p) const { return sign(r - (c - p).mag()) >= 0; }
	ld area(const ld& lo, const ld& hi) const { return (hi - lo) * r * r * .5; }
	ld rad(const Pos& p) const { return (p - c).rad(); }
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld green(const ld& lo, const ld& hi) const {
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
};
ld intersection(const Seg& s1, const Seg& s2) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	return a1;
	//if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	//return -1;
}
//Vld circle_line_intersections(const Pos& s, const Pos& e, const Circle& q, const bool& f = 0) {
//	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
//	Pos vec = e - s;
//	Pos OM = s - q.c;
//	ld a = vec.Euc();
//	ld b = vec * OM;
//	ld c = OM.Euc() - q.r * q.r;
//	ld J = b * b - a * c;
//	if (J < -TOL) return {};
//	ld det = sqrt(std::max((ld)0, J));
//	ld lo = (-b - det) / a;
//	ld hi = (-b + det) / a;
//	Vld ret;
//	if (f) {
//		if (0 < hi && hi < 1) ret.push_back(hi);
//		if (zero(det)) return ret;
//		if (0 < lo && lo < 1) ret.push_back(lo);
//	}
//	else {
//		auto the = [&](ld rt) { return q.rad(s + (e - s) * rt); };
//		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
//		if (zero(det)) return ret;
//		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
//	}
//	return ret;
//}
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
		//if (-TOL < hi && hi < 1 + TOL) ret.push_back(hi);
		//if (zero(det)) return ret;
		//if (-TOL < lo && lo < 1 + TOL) ret.push_back(lo);
		ret.push_back(hi); ret.push_back(lo);
	}
	else {
		auto the = [&](const ld& rt) { return norm(q.rad(s + (e - s) * rt)); };
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
		if (zero(det)) return ret;
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
	}
	return ret;
}
Vld intersection(const Circle& a, const Circle& b) {
	Pos ca = a.c, cb = b.c;
	Pos vec = cb - ca;
	ld ra = a.r, rb = b.r;
	ld distance = vec.mag();
	ld rd = vec.rad();
	if (vec.Euc() > sq(ra + rb) + TOL) return {};
	if (vec.Euc() < sq(ra - rb) - TOL) return {};
	ld X = (ra * ra - rb * rb + vec.Euc()) / (2 * distance * ra);
	if (X < -1) X = -1;
	if (X > 1) X = 1;
	ld h = acos(X);
	Vld ret = {};
	ret.push_back(norm(rd + h));
	if (zero(h)) return ret;
	ret.push_back(norm(rd - h));
	return ret;
}
struct Arc {
	ld lo, hi;
	Arc(ld l_ = 0, ld h_ = 0) : lo(l_), hi(h_) {}
	bool operator < (const Arc& a) const { return zero(lo - a.lo) ? hi < a.hi : lo < a.lo; }
	inline friend std::istream& operator >> (std::istream& is, Arc& a) { is >> a.lo >> a.hi; return is; }
	inline friend std::ostream& operator << (std::ostream& os, const Arc& a) { os << a.lo << " " << a.hi; return os; }
} arc[LEN]; bool F[LEN];
typedef std::vector<Arc> Arcs;
Circle C[LEN];
Polygon B[LEN];
Seg S[LEN];
int shadow(const Circle& c0, const Circle& c1, Arc& a1, Arc& a2) {
	Vld inxs = intersection(c0, c1);
	if (inxs.size() < 2) return 0;
	ld lo = inxs[0], hi = inxs[1];
	if (lo < hi) { a1 = Arc(lo, hi); return 1; }
	a1 = Arc(lo, 2 * PI);
	a2 = Arc(0, hi);
	return 2;
}
//int shadow(const Circle& c, const Polygon& b, Arc& a1, Arc& a2) {
//	ld lo = -1, hi = -1;
//	assert(b.size() == 4);
//	int f = 0;
//	for (int i = 0, j; i < 4; i++) {
//		j = (i + 1) % 4;
//		const Pos& p0 = b[i], & p1 = b[j];
//		if (sign(cross(p0, p1, c.c) / (p0 - p1).mag() - c.r) >= 0) f++;
//	}
//	if (f >= 4) { a1 = Arc(0, 2 * PI); return 1; }
//	bool f0 = 0, f1 = 0;
//	for (int i = 0, j; i < 4; i++) {
//		j = (i + 1) % 4;
//		const Pos& p0 = b[j], & p1 = b[i];
//		Vld inxs = circle_line_intersections(c, p0, p1, CIRCLE);
//		int sz = inxs.size();
//		if (sz == 0) continue;
//		if (sz == 2) { f0 = f1 = 1; lo = inxs[0]; hi = inxs[1]; break; }
//		if (c >= p0) f0 = 1, hi = inxs[0];
//		else if (c >= p1) f1 = 1, lo = inxs[1];
//	}
//	if (!f0 || !f1) return 0;
//	assert(-TOL < lo && lo < 2 * PI + TOL);
//	assert(-TOL < hi && hi < 2 * PI + TOL);
//	if (lo < hi) { a1 = Arc(lo, hi); return 1; }
//	a1 = Arc(lo, 2 * PI); a2 = Arc(0, hi);
//	return 2;
//}
int shadow(const Circle& c, const Polygon& b, Arcs& va) {
	assert(b.size() == 4);
	int f = 0;
	for (int i = 0, j; i < 4; i++) {
		j = (i + 1) % 4;
		const Pos& p0 = b[i], & p1 = b[j];
		if (sign(cross(p0, p1, c.c) / (p0 - p1).mag() - c.r) >= 0) f++;
	}
	if (f >= 4) { va = { Arc(0, 2 * PI) }; return 1; }
	Vld vx = { 0, 2 * PI };
	for (int i = 0, j; i < 4; i++) {
		j = (i + 1) % 4;
		const Pos& p0 = b[j], & p1 = b[i];
		Vld inxs = circle_line_intersections(c, p0, p1, CIRCLE);
		for (const ld& x : inxs) vx.push_back(x);
	}
	std::sort(vx.begin(), vx.end());
	vx.erase(unique(vx.begin(), vx.end()), vx.end());
	int sz = vx.size();
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		ld lo = vx[i], hi = vx[j];
		ld x = (vx[i] + vx[j]) * .5;
		Pos mid = c.p(x);
		if (inner_check(b, mid)) va.push_back(Arc(lo, hi));
	}
	return 1;
}
int shadow(const Seg& s, const Circle& c, Pos& p) {
	Vld inxs = circle_line_intersections(c, s.s, s.e, LINE);
	int sz = inxs.size();
	if (sz < 2) return 0;
	ld lo = inxs[0], hi = inxs[1];
	lo = fit(lo); hi = fit(hi);
	if (zero(lo - hi)) return 0;
	p = Pos(lo, hi);
	return 1;
}
int shadow(const Seg& s, const Polygon& b, Pos& p) {
	int sz = b.size();
	assert(sz == 4);
	ld lo = -1, hi = -1;
	int f = 0;
	for (int i = 0; i < sz; i++) if (ccw(s.s, s.e, b[i]) >= 0) f++;
	if (f >= sz) return 0;
	bool f0 = 0, f1 = 0;
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		const Pos& p0 = b[i], & p1 = b[j];
		if (!ccw(s.s, s.e, p0, p1)) continue;
		if (ccw(s.s, s.e, p0) * ccw(s.s, s.e, p1) <= 0) {
			ld ix = intersection(s, Seg(p0, p1));
			if (ccw(s.s, s.e, p0, p1) > 0) f0 = 1, hi = ix;
			else f1 = 1, lo = ix;
		}
	}
	if (!f0 || !f1) return 0;
	lo = fit(lo); hi = fit(hi);
	if (zero(lo - hi)) return 0;
	p = Pos(lo, hi);
	return 1;
}
bool remain(const Polygon& P, const ld& m, Pos& ret) {
	auto box = [&](const int& i, const int& j) -> Polygon {
		Pos v = ~(P[j] - P[i]).unit() * m;
		return { P[j] + v, P[i] + v, P[i] - v, P[j] - v };
	};
	int sz = P.size();
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		C[i] = Circle(P[i], m);
		B[i] = box(i, j);
		S[i] = Seg(B[i][0], B[i][1]);
	}
	for (int i = 0; i < sz; i++) {
		if (!F[i]) continue;
		Arcs VA;
		Arc a = arc[i];
		if (a.lo < a.hi) VA = { a };
		else VA = { Arc(a.lo, 2 * PI), Arc(0, a.hi) };
		for (int j = 0; j < sz; j++) {
			if (i == j) continue;
			Arc a1, a2;
			int f;
			f = shadow(C[i], C[j], a1, a2);
			if (f) VA.push_back(a1);
			if (f == 2) VA.push_back(a2);
			if (i == (j + 1) % sz) continue;
			Arcs va;
			f = shadow(C[i], B[j], va);
			for (const Arc& a : va) VA.push_back(a);
			//f = shadow(C[i], B[j], a1, a2);
			//if (f) VA.push_back(a1);
			//if (f == 2) VA.push_back(a2);
		}
		std::sort(VA.begin(), VA.end());
		VA.push_back(Arc(2 * PI, 2 * PI));
		ld hi = 0;
		for (const Arc& a : VA) {
			if (hi < a.lo) { ret = C[i].p(hi); return 1; }
			else hi = std::max(hi, a.hi);
		}
	}
	for (int i = 0; i < sz; i++) {
		Polygon VP;
		for (int j = 0; j < sz; j++) {
			if (i == j) continue;
			Pos p;
			int f;
			f = shadow(S[i], B[j], p);
			if (f) VP.push_back(p);
			if ((i + 1) % sz == j) continue;
			f = shadow(S[i], C[j], p);
			if (f) VP.push_back(p);
		}
		std::sort(VP.begin(), VP.end());
		VP.push_back(Pos(1, 1));
		ld hi = 0;
		for (const Pos& p : VP) {
			if (hi < p.LO) { ret = S[i].p(hi); return 1; }
			else hi = std::max(hi, p.HI);
		}
	}
	return 0;
}
ld bi_search(const Polygon& P, Pos& ret) {
	int sz = P.size();
	ld s = TOL, e = -1;
	for (int i = 0; i < sz; i++) {
		for (int j = i + 1; j < N; j++) {
			ld d = (P[i] - P[j]).mag();
			e = std::max(e, d);
		}
	}
	e *= .5; e += TOL;
	std::cout << "DEBUG::\n";
	std::cout << "	s:: " << s << " e:: " << e << "\n";
	std::cout << "DEBUG::\n";
	int c = 50; while (c--) {
		ld m = (s + e) * .5;
		if (remain(P, m, ret)) s = m;
		else e = m;
	}
	return (s + e) * .5 / V;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	std::cin >> V >> N;
	Polygon P(N); for (Pos& p : P) std::cin >> p; norm(P);
	for (int i = 0, i0, i1, i2; i < N; i++) {
		i0 = (i - 1 + N) % N;
		i1 = i;
		i2 = (i + 1) % N;
		if (ccw(P[i0], P[i1], P[i2]) >= 0) { F[i] = 0; continue; }
		Pos v;
		v = ~(P[i1] - P[i0]);
		arc[i].lo = v.rad();
		v = ~(P[i2] - P[i1]);
		arc[i].hi = v.rad();
		F[i] = 1;
	}
	Pos ret;
	std::cout << bi_search(P, ret) << "\n";
	std::cout << ret << "\n";
	return;
}
int main() { solve(); return 0; }//boj35239
