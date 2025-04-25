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
const ld TOL = 1e-7;
const ld PI = acos(-1);
const int LEN = 1e3;
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
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d2, d3)) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d2, d3)) > 0; }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
ld dist(const Pos& d1, const Pos& d2, const Pos& t, bool f = 0) {
	if (!f) return cross(d1, d2, t) / (d1 - d2).mag();
	if (sign(projection(d1, d2, d2, t)) <= 0 &&
		sign(projection(d2, d1, d1, t)) <= 0)
		return std::abs(cross(d1, d2, t)) / (d1 - d2).mag();
	return std::min((d1 - t).mag(), (d2 - t).mag());
}
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
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
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	std::sort(C.begin(), C.end());
	if (C.size() <= 2) { for (const Pos& p : C) H.push_back(p); }
	else {
		for (int i = 0; i < C.size(); i++) {
			while (H.size() > 1 && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
		int s = H.size() + 1;
		for (int i = C.size() - 1; i >= 0; i--) {
			while (H.size() > s && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
	}
	return H;
}
Polygon convex_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
	Polygon qs;
	int n = ps.size();
	for (int i = 0; i < n; i++) {
		Pos p1 = ps[i], p2 = ps[(i + 1) % n];
		int d1 = ccw(b1, b2, p1), d2 = ccw(b1, b2, p2);
		if (d1 >= 0) qs.push_back(p1);
		if (d1 * d2 < 0) qs.push_back(intersection(p1, p2, b1, b2));
	}
	return qs;
}
Polygon sutherland_hodgman(const Polygon& C, const Polygon& clip) {
	int sz = clip.size();
	std::vector<Pos> ret = C;
	for (int i = 0; i < sz; i++) {
		Pos b1 = clip[i], b2 = clip[(i + 1) % sz];
		ret = convex_cut(ret, b1, b2);
	}
	return ret;
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
Segs half_plane_intersection(Segs& HP, const bool& srt = 1) {
	auto check = [&](Seg& u, Seg& v, Seg& w) -> bool {
		return w.inner(intersection_(u, v));
		};
	if (srt) std::sort(HP.begin(), HP.end());
	std::deque<Seg> dq;
	int sz = HP.size();
	for (int i = 0; i < sz; ++i) {
		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
		dq.push_back(HP[i]);
	}
	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
	sz = dq.size();
	if (sz < 3) return {};
	std::vector<Seg> HPI;
	for (int i = 0; i < sz; ++i) HPI.push_back(dq[i]);
	return HPI;
}
Segs half_plane_intersection(const Segs& P, const Segs& Q) {
	Segs HP(P.size() + Q.size());
	std::merge(P.begin(), P.end(), Q.begin(), Q.end(), HP.begin());
	return half_plane_intersection(HP, 0);
}

struct Circle {
	Pos c;
	int r;
	Circle(Pos c_ = Pos(), int r_ = 0) : c(c_), r(r_) {}
	bool operator == (const Circle& q) const { return c == q.c && r == q.r; }
	bool operator != (const Circle& q) const { return !(*this == q); }
	bool operator < (const Circle& q) const { return c == q.c ? r < q.r : c < q.c; }
	//bool operator < (const Circle& q) const { return r < q.r && (c - q.c).mag() + r < q.r + TOL; }
	bool outside(const Circle& q) const { return sign((c - q.c).Euc() - sq((ll)r + q.r)) >= 0; }
	Circle operator + (const Circle& q) const { return { c + q.c, r + q.r }; }
	Circle operator - (const Circle& q) const { return { c - q.c, r - q.r }; }
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld rad(const Pos& p) const { return (p - c).rad(); }
	ld area(const ld& lo, const ld& hi) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
	ld H(const ld& th) const { return sin(th) * c.x + cos(th) * c.y + r; }//coord trans | check right
	//bool operator < (const Pos& p) const { return r < (c - p).mag(); }
	bool operator < (const Pos& p) const { return sign(r - (c - p).mag()) < 0; }
	bool operator > (const Pos& p) const { return r > (c - p).mag(); }
	bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
	friend std::istream& operator >> (std::istream& is, Circle& c) { is >> c.c >> c.r; return is; }
	friend std::ostream& operator << (std::ostream& os, const Circle& c) { os << c.c << " " << c.r; return os; }
} INVAL = { { 0, 0 }, -1 };
bool cmpr(const Circle& p, const Circle& q) { return p.r > q.r; }//sort descending order
Vld intersections(const Circle& a, const Circle& b) {
	Pos ca = a.c, cb = b.c;
	Pos vec = cb - ca;
	ll ra = a.r, rb = b.r;
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
Circle enclose_circle(const Pos& u, const Pos& v) {
	Pos c = (u + v) * .5;
	return Circle(c, (c - u).mag());
}
//Circle enclose_circle(const Pos& u, const Pos& v, const Pos& w) {
//	Line l1 = rotate90(L(u, v), (u + v) * .5);
//	Line l2 = rotate90(L(v, w), (v + w) * .5);
//	if (zero(l1 / l2)) return { { 0, 0 }, -1 };
//	Pos c = intersection(l1, l2);
//	ld r = (c - u).mag();
//	return Circle(c, r);
//}
Circle enclose_circle(const Pos& u, const Pos& v, const Pos& w) {
	if (!ccw(u, v, w)) return INVAL;
	Pos m1 = (u + v) * .5, v1 = ~(v - u);
	Pos m2 = (u + w) * .5, v2 = ~(w - u);
	Pos c = intersection(m1, m1 + v1, m2, m2 + v2);
	return Circle(c, (u - c).mag());
}
Circle enclose_circle(std::vector<Pos> R) {
	if (R.size() == 0) return Circle(O, -1);
	else if (R.size() == 1) return Circle(R[0], 0);
	else if (R.size() == 2) return enclose_circle(R[0], R[1]);
	else return enclose_circle(R[0], R[1], R[2]);
}
Circle minimum_enclose_circle(std::vector<Pos> P) {
	shuffle(P.begin(), P.end(), std::mt19937(0x14004));
	Circle mec = INVAL;
	int sz = P.size();
	for (int i = 0; i < sz; i++) {
		if (mec.r < -1 || mec < P[i]) {
			mec = Circle(P[i], 0);
			for (int j = 0; j <= i; j++) {
				if (mec < P[j]) {
					Circle ans = enclose_circle(P[i], P[j]);
					if (zero(mec.r)) { mec = ans; continue; }
					Circle l = INVAL, r = INVAL;
					//Pos vec = P[j] - P[i];
					for (int k = 0; k <= j; k++) {
						if (ans < P[k]) {
							//ld CCW = vec / (P[k] - P[j]);
							ld CCW = cross(P[i], P[j], P[k]);
							Circle c = enclose_circle(P[i], P[j], P[k]);
							if (c.r < 0) continue;
							//else if (CCW > 0 && (l.r < 0 || (vec / (c.c - P[i])) > (vec / (l.c - P[i])))) l = c;
							//else if (CCW < 0 && (r.r < 0 || (vec / (c.c - P[i])) < (vec / (r.c - P[i])))) r = c;
							else if (CCW > 0 && (l.r < 0 || cross(P[i], P[j], c.c) > cross(P[i], P[j], l.c))) l = c;
							else if (CCW < 0 && (r.r < 0 || cross(P[i], P[j], c.c) < cross(P[i], P[j], r.c))) r = c;
						}
					}
					if (l.r < 0 && r.r < 0) mec = ans;
					else if (l.r < 0) mec = r;
					else if (r.r < 0) mec = l;
					else mec = l.r < r.r ? l : r;
				}
			}
		}
	}
	return mec;
}

Vld circle_line_intersections(const Seg& l, const Circle& q, const int& t = LINE) {
	//https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
	Pos s = l.s, e = l.e;
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
	if (t == LINE) {
		if (0 < lo && lo < 1) ret.push_back(lo);
		if (zero(det)) return ret;
		if (0 < hi && hi < 1) ret.push_back(hi);
	}
	else {//circle
		auto the = [&](ld rt) { return q.rad(s + (e - s) * rt); };
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
		if (zero(det)) return ret;
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
	}
	return ret;
}
Vld circle_line_intersections(const Pos& s, const Pos& e, const Circle& q, const bool& f = 0) {
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
	if (f) {
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
Polygon circle_line_intersection(const Pos& o, const ld& r, const Pos& p1, const Pos& p2) {
	ld d = dist(p1, p2, o);
	if (std::abs(d) > r) return {};
	Pos vec = p2 - p1;
	Pos m = intersection(p1, p2, o, o + ~vec);
	ld distance = vec.mag();
	ld ratio = sqrt(r * r - d * d);
	Pos m1 = m - vec * ratio / distance;
	Pos m2 = m + vec * ratio / distance;
	if (dot(p1, p2, m1, m2) < 0) std::swap(m1, m2);
	return { m1, m2 };//p1->p2
}
ld circle_cut(const Circle& c, const Pos& p1, const Pos& p2) {
	Pos v1 = p1 - c.c, v2 = p2 - c.c;
	ld r = c.r;
	std::vector<Pos> inx = circle_line_intersection(O, r, v1, v2);
	if (inx.empty()) return r * r * rad(v1, v2) * .5;
	Pos m1 = inx[0], m2 = inx[1];
	bool d1 = dot(m1, v1, m2) > -TOL, d2 = dot(m1, v2, m2) > -TOL;
	if (d1 && d2) return (v1 / v2) * .5;
	else if (d1) return (v1 / m2 + r * r * rad(m2, v2)) * .5;
	else if (d2) return (r * r * rad(v1, m1) + m1 / v2) * .5;
	else if (dot(v1, m1, v2) > 0 && dot(v1, m2, v2) > 0)
		return (r * r * (rad(v1, m1) + rad(m2, v2)) + m1 / m2) * .5;
	else return (r * r * rad(v1, v2)) * .5;
}

struct Arc {
	ld lo, hi;// [lo, hi] - radian range of arc
	Circle c;// c.r - radius of arc
	Arc(ld l_ = 0, ld h_ = 0) : lo(l_), hi(h_) {}
	bool operator < (const Arc& a) const { return zero(lo - a.lo) ? hi < a.hi : lo < a.lo; }
	inline friend std::istream& operator >> (std::istream& is, Arc& a) { is >> a.lo >> a.hi; return is; }
	inline friend std::ostream& operator << (std::ostream& os, const Arc& a) { os << a.lo << " " << a.hi; return os; }
};

struct Triangle {
	Pos a, b, c;
	inline Triangle(Pos p = Pos(), Pos q = Pos(), Pos r = Pos()) {
		if ((q - p) / (r - q) < 0) std::swap(q, r);
		a = p; b = q; c = r;
	}
	inline int inner_check(const Pos& p, const Pos& v) const {
		ld f1 = (b - a) / (p - a);
		ld f2 = (c - b) / (p - b);
		ld f3 = (a - c) / (p - c);
		if (sign(f1) < 0 || sign(f2) < 0 || sign(f3) < 0) return 0;
		if (zero(f1)) return (b - a) / v > 0 ? 2 : 0;//on_seg && centripetal
		if (zero(f2)) return (c - b) / v > 0 ? 2 : 0;
		if (zero(f3)) return (a - c) / v > 0 ? 2 : 0;
		return 1;
	}
};

struct Vec {
	ld vy, vx;
	Vec(ld Y = 0, ld X = 0) : vy(Y), vx(X) {}
	bool operator == (const Vec& v) const { return (zero(vy - v.vy) && zero(vx - v.vx)); }
	bool operator < (const Vec& v) const { return zero(vy - v.vy) ? vx < v.vx : vy < v.vy; }
	ld operator * (const Vec& v) const { return vy * v.vy + vx * v.vx; }
	ld operator / (const Vec& v) const { return vy * v.vx - vx * v.vy; }
	Vec operator ~ () const { return { -vx, vy }; }
	Vec& operator *= (const ld& n) { vy *= n; vx *= n; return *this; }
	Vec& operator /= (const ld& n) { vy /= n; vx /= n; return *this; }
	ld mag() const { return hypot(vy, vx); }
}; const Vec Zero = { 0, 0 };
struct Line {//ax + by = c
	Vec s;
	ld c;
	Line(Vec V = Vec(0, 0), ld C = 0) : s(V), c(C) {}
	bool operator < (const Line& l) const {
		bool f1 = Zero < s;
		bool f2 = Zero < l.s;
		if (f1 != f2) return f1;
		ld CCW = s / l.s;
		return zero(CCW) ? c * hypot(l.s.vy, l.s.vx) < l.c * hypot(s.vy, s.vx) : CCW > 0;
	}
	ld operator * (const Line& l) const { return s * l.s; }
	ld operator / (const Line& l) const { return s / l.s; }
	Line operator + (const ld& n) const { return Line(s, c + hypot(s.vy, s.vx) * n); }
	Line operator - (const ld& n) const { return Line(s, c - hypot(s.vy, s.vx) * n); }
	Line operator * (const ld& n) const { return Line({ s.vy * n, s.vx * n }, c * n); }
	Line& operator += (const ld& n) { c += hypot(s.vy, s.vx) * n; return *this; }
	Line& operator -= (const ld& n) { c -= hypot(s.vy, s.vx) * n; return *this; }
	Line& operator *= (const ld& n) { s *= n, c *= n; return *this; }
	ld dist(const Pos& p) const { return s.vy * p.x + s.vx * p.y; }
	ld above(const Pos& p) const { return s.vy * p.x + s.vx * p.y - c; }
	ld mag() const { return s.mag(); }
	friend inline ld rad(const Line& b, const Line& l) { return atan2(b / l, b * l); }
	friend std::ostream& operator << (std::ostream& os, const Line& l) { os << l.s.vy << " " << l.s.vx << " " << l.c; return os; }
};
const Line Xaxis = { { 0, -1 }, 0 };
const Line Yaxis = { { 1, 0 }, 0 };
Line L(const Pos& s, const Pos& e) {
	ld dy, dx, c;
	dy = e.y - s.y;
	dx = s.x - e.x;
	c = dy * s.x + dx * s.y;
	return Line(Vec(dy, dx), c);
}
Line L(const Vec& s, const Pos& p) {
	ld c = s.vy * p.x + s.vx * p.y;
	return Line(s, c);
}
Line rotate(const Line& l, const Pos& p, ld the) {
	Vec s = l.s;
	ld x = -s.vx, y = s.vy;
	ld vx = -(x * cos(the) - y * sin(the));
	ld vy = x * sin(the) + y * cos(the);
	ld c = vy * p.x + vx * p.y;
	return Line(Vec(vy, vx), c);
}
Line rot90(const Line& l, const Pos& p) {
	Vec s = ~l.s;
	ld c = s.vy * p.x + s.vx * p.y;
	return Line(s, c);
}
Pos intersection(const Line& l1, const Line& l2) {
	Vec v1 = l1.s, v2 = l2.s;
	ld det = v1 / v2;
	return {
		(l1.c * v2.vx - l2.c * v1.vx) / det,
		(l2.c * v1.vy - l1.c * v2.vy) / det,
	};
}
ld rad(const Line& b, const Line& l) { return atan2(b / l, b * l); }
bool half_plane_intersection(std::vector<Line>& HP, std::vector<Pos>& hull) {
	auto cw = [&](const Line& l1, const Line& l2, const Line& target) -> bool {
		if (l1.s / l2.s < TOL) return 0;
		Pos p = intersection(l1, l2);
		//return target.s.vy * p.x + target.s.vx * p.y > target.c - TOL;
		return target.above(p) > -TOL;
		};
	std::deque<Line> dq;
	std::sort(HP.begin(), HP.end());
	for (const Line& l : HP) {
		if (!dq.empty() && zero(dq.back() / l)) continue;
		while (dq.size() >= 2 && cw(dq[dq.size() - 2], dq.back(), l)) dq.pop_back();
		while (dq.size() >= 2 && cw(l, dq.front(), dq[1])) dq.pop_front();
		dq.push_back(l);
	}
	while (dq.size() >= 3 && cw(dq[dq.size() - 2], dq.back(), dq.front())) dq.pop_back();
	while (dq.size() >= 3 && cw(dq.back(), dq.front(), dq[1])) dq.pop_front();
	for (int i = 0; i < (int)dq.size(); i++) {
		Line cur = dq[i], nxt = dq[(i + 1) % (int)dq.size()];
		if (cur / nxt < TOL) {
			hull.clear();
			return 0;
		}
		hull.push_back(intersection(cur, nxt));
	}
	return 1;
}

struct Linear {//ps[0] -> ps[1] :: refer to bulijiojiodibuliduo
	Pos ps[2];
	Pos dir_;
	Pos& operator[](int i) { return ps[i]; }
	Pos dir() const { return dir_; }
	Linear(Pos a = Pos(0, 0), Pos b = Pos(0, 0)) {
		ps[0] = a;
		ps[1] = b;
		dir_ = (ps[1] - ps[0]).unit();
	}
	bool include(const Pos& p) const { return sign(dir_ / (p - ps[0])) > 0; }
	Linear push() const {//push eps outward
		const double eps = 1e-8;
		Pos delta = ~(ps[1] - ps[0]).unit() * eps;
		return Linear(ps[0] + delta, ps[1] + delta);
	}
	Linear operator + (const double eps) const {//push eps outward
		Pos delta = ~(ps[1] - ps[0]).unit() * eps;
		return Linear(ps[0] + delta, ps[1] + delta);
	}
	Linear operator - (const double eps) const {//pull eps inward
		Pos delta = ~(ps[1] - ps[0]).unit() * eps;
		return Linear(ps[0] - delta, ps[1] - delta);
	}
	friend bool parallel(const Linear& l0, const Linear& l1) { return zero(l0.dir() / l1.dir()); }
	friend bool same_dir(const Linear& l0, const Linear& l1) { return parallel(l0, l1) && l0.dir() * l1.dir() > 0; }
	bool operator < (const Linear& l0) const {
		if (same_dir(*this, l0)) return l0.include(ps[0]);
		else return cmpq(this->dir(), l0.dir());
	}
};
Pos intersection(Linear& l1, Linear& l2) { return intersection(l1[0], l1[1], l2[0], l2[1]); }
std::vector<Pos> half_plane_intersection(std::vector<Linear>& HP) {//refer to bulijiojiodibuliduo
	auto check = [&](Linear& u, Linear& v, Linear& w) -> bool {
		return w.include(intersection(u, v));
		};
	std::sort(HP.begin(), HP.end());
	std::deque<Linear> dq;
	int sz = HP.size();
	for (int i = 0; i < sz; ++i) {
		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
		dq.push_back(HP[i]);
	}
	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
	sz = dq.size();
	if (sz < 3) return {};
	std::vector<Pos> HPI;
	for (int i = 0; i < sz; ++i) HPI.push_back(intersection(dq[i], dq[(i + 1) % sz]));
	return HPI;
}



//3D============================================================================//

struct Pos3D {
	ld x, y, z;
	int r, i;
	Pos3D(ld x_ = 0, ld y_ = 0, ld z_ = 0) : x(x_), y(y_), z(z_) {}
	bool operator == (const Pos3D& p) const { return zero(x - p.x) && zero(y - p.y) && zero(z - p.z); }
	bool operator != (const Pos3D& p) const { return !zero(x - p.x) || !zero(y - p.y) || !zero(z - p.z); }
	bool operator < (const Pos3D& p) const { return zero(x - p.x) ? zero(y - p.y) ? z < p.z : y < p.y : x < p.x; }
	ld operator * (const Pos3D& p) const { return x * p.x + y * p.y + z * p.z; }
	Pos3D operator / (const Pos3D& p) const { return { y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x }; }
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	Pos3D operator * (const ld& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ld& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ld& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ld& n) { x /= n; y /= n; z /= n; return *this; }
	ld Euc() const { return x * x + y * y + z * z; }
	ld mag() const { return sqrtl(Euc()); }
	ld lon() const { return atan2(y, x); }
	ld lat() const { return atan2(z, sqrtl(x * x + y * y)); }
	Pos3D operator - () const { return { -x, -y, -z }; }
	Pos3D unit() const { return *this / mag(); }
	Pos3D norm(const Pos3D& p) const { return (*this / p).unit(); }
	Pos3D rodrigues_rotate(const ld& th, const Pos3D& axis) const {
		ld s = sin(th), c = cos(th); Pos3D n = axis.unit();
		return n * (*this * n) * (1 - c) + (*this * c) + (n / *this) * s;
	}
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z; return os; }
};
typedef std::vector<Pos3D> Polyhedron;
const Pos3D O3D = { 0, 0, 0 };
const Pos3D X_axis = { 1, 0, 0 };
const Pos3D Y_axis = { 0, 1, 0 };
const Pos3D Z_axis = { 0, 0, 1 };
const Pos3D MAXP3D = { INF, INF, INF };
std::vector<Pos3D> pos;
Pos3D S2C(const ld& lon, const ld& lat) {//Spherical to Cartesian
	ld phi = lon * PI / 180;
	ld the = lat * PI / 180;
	return Pos3D(cos(phi) * cos(the), sin(phi) * cos(the), sin(the));
}
Pos3D point(const Pos3D Xaxis, const Pos3D Yaxis, const ld& th) { return Xaxis * cos(th) + Yaxis * sin(th); }
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) { return norm(atan2(Yaxis * p, Xaxis * p)); }
int above(const Plane& S, const Pos3D& p) { return sign(p * S.norm() + S.d); }
//ld randTOL() {
//	ld rand01 = rand() / (ld)RAND_MAX;
//	ld err = (rand01 - .5) * TOL;
//	return err;
//}
//Pos3D add_noise(const Pos3D& p) {//refer to BIGINTEGER
//	return p + Pos3D(randTOL(), randTOL(), randTOL());
//}
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
ld dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
ld area(const std::vector<Pos3D>& H, const Pos3D& norm) {
	if (H.size() < 3) return 0;
	ld ret = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos3D& cur = H[i], nxt = H[(i + 1) % sz];
		ret += cross(H[0], cur, nxt) * norm / norm.mag();
	}
	return std::abs(ret * .5);
}
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) > 0; }
//std::vector<Pos3D> graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
ld graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
	//if (C.size() < 3) {
	//	std::sort(C.begin(), C.end());
	//	return C;
	// }
	if (C.size() < 3) return 0;
	std::vector<Pos3D> H;
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pos3D& p, const Pos3D& q) -> bool {
		ld ret = ccw(C[0], p, q, norm);
		if (zero(ret)) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i], norm) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	//return H;
	return area(H, norm);
}
bool inner_check(const Pos3D& d1, const Pos3D& d2, const Pos3D& t) {
	Pos3D nrm = cross(O3D, d1, d2);
	Pos3D p1 = d1, p2 = d2;
	if (ccw(O3D, p1, p2, nrm) < 0) std::swap(p1, p2);
	return ccw(O3D, p1, t, nrm) >= 0 && ccw(O3D, p2, t, nrm) <= 0;
}
int inner_check(std::vector<Pos3D>& H, const Pos3D& p) {//for convex hull
	int sz = H.size();
	if (sz <= 1) return -1;
	if (sz == 2) {
		if (on_seg_strong(H[0], H[1], p)) return 0;
		else return -1;
	}
	Pos3D torque0 = cross(H[0], H[1], p);
	for (int i = 1; i < sz; i++) {
		Pos3D cur = H[i], nxt = H[(i + 1) % sz];
		Pos3D torqueI = cross(cur, nxt, p);
		if (zero(torqueI.mag())) {
			if (on_seg_strong(cur, nxt, p)) return 0;
			else return -1;
		}
		if (torque0 * torqueI < 0) return -1;
	}
	return 1;
}

struct Line3D {
	Pos3D dir, p0;
	Line3D(Pos3D dir_ = Pos3D(0, 0, 0), Pos3D p0_ = Pos3D(0, 0, 0)) : dir(dir_), p0(p0_) {}
};
Line3D L(const Pos3D& p1, const Pos3D& p2) { return { p2 - p1, p1 }; }
Line3D line(const Pos3D& p1, const Pos3D& p2) { return { p2 - p1, p1 }; }

struct Plane {
	ld a, b, c, d;
	Plane(ld a_ = 0, ld b_ = 0, ld c_ = 0, ld d_ = 0) : a(a_), b(b_), c(c_), d(d_) {}
	Pos3D norm() const { return Pos3D(a, b, c); };
	Plane operator + (const ld& n) const { return { a, b, c, d + n }; }
	Plane operator - (const ld& n) const { return { a, b, c, d - n }; }
	Plane& operator += (const ld& n) { d += n; return *this; }
	Plane& operator -= (const ld& n) { d -= n; return *this; }
	friend std::istream& operator >> (std::istream& is, Plane& f) { is >> f.a >> f.b >> f.c >> f.d; return is; }
	friend std::ostream& operator << (std::ostream& os, const Plane& f) { os << f.a << " " << f.b << " " << f.c << " " << f.d; return os; }
} knife;
Plane plane(const Pos3D& p, const ld& n) { return Plane(p.x, p.y, p.z, n); }
ld dist(const Plane& s, const Pos3D& p) { return (s.norm() * p + s.d) / s.norm().mag(); }
Pos3D intersection(const Plane& S, const Line3D& l) {
	ld det = S.norm() * l.dir;
	if (zero(det)) return { INF, INF, INF };
	ld t = -((S.norm() * l.p0 + S.d) / det);
	return l.p0 + (l.dir * t);
}
Pos3D intersection(const Plane& S, const Pos3D& p1, const Pos3D& p2, const bool& f = 0) {
	Line3D l = line(p1, p2);
	Pos3D inx = intersection(S, l);
	if (f && !on_seg_strong(p1, p2, inx)) return { INF, INF, INF };
	return inx;
}
bool circle_intersection(const Pos3D& a, const Pos3D& b, const ld& th, std::vector<Pos3D>& inxs) {
	inxs.clear();
	Pos3D mid = (a + b) * .5;
	if (zero(mid.mag())) return 0;
	ld x = cos(th) / mid.mag();
	if (x < -1 || 1 < x) return 0;
	Pos3D w = mid.unit() * x;
	ld ratio = sqrtl(1 - x * x);
	Pos3D h = (mid.unit() / (b - a).unit()) * ratio;
	inxs.push_back(w + h);
	inxs.push_back(w - h);
	return 1;
}
//bool plane_circle_intersection(const Pos3D& a, const Pos3D& perp, const ld& th, std::vector<Pos3D>& inxs) {
//	inxs.clear();
//	Pos3D vec = a - (perp * (perp * a));
//	if (zero(vec.mag())) return 0;
//	ld x = cos(th) / vec.mag();
//	if (x < -1 || 1 < x) return 0;
//	Pos3D w = vec.unit() * x;
//	ld ratio = sqrtl(1 - x * x);
//	Pos3D h = (vec.unit() / perp) * ratio;
//	inxs.push_back(w + h);
//	inxs.push_back(w - h);
//	return 1;
//}
int intersection(const Plane& p1, const Plane& p2, Line3D& l) {
	Pos3D n1 = p1.norm();
	Pos3D n2 = p2.norm();
	Pos3D dir = n2 / n1;
	if (zero(dir.mag())) {
		ld f = n1 * n2;
		ld d1 = dist(p1, O3D);
		ld d2 = dist(p2, O3D);
		if (sign(f) > 0) return sign(d2 - d1) >= 0 ? 0 : -1;
		else return sign(d2 + d1) >= 0 ? 0 : -2;
	}
	dir = dir.unit();
	Pos3D q1 = intersection(p1, Line3D(n1, O3D));
	Pos3D v1 = n1 / dir;
	Pos3D p0 = intersection(p2, Line3D(v1, q1));
	l = Line3D(dir, p0);
	return 1;
}
Pos3D s2c(const ld& phi, const ld& psi) {//Spherical to Cartesian
	ld lat = phi * PI / 180, lon = psi * PI / 180;
	return Pos3D(cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat));
}
bool plane_circle_intersection(const Pos3D& perp, const Pos3D& a, Polyhedron& inxs, const ld& R = 0) {
	inxs.clear();
	Pos3D vec = a - (perp * (perp * a));
	if (zero(vec.mag())) return 0;
	ld th = (ld)a.r / R;
	ld x = cos(th) / vec.mag();
	if (x < -1 || 1 < x) return 0;
	Pos3D w = vec.unit() * x;
	ld ratio = sqrtl(1 - x * x);
	Pos3D h = (vec.unit() / perp) * ratio;
	inxs.push_back(w + h);
	if (!zero(ratio)) inxs.push_back(w - h);
	return 1;
}
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) {
	ld X = Xaxis * p, Y = Yaxis * p;
	ld th = atan2(Y, X);
	return th;
}
ld angle(const Pos3D& a, const Pos3D& b) {
	Pos3D perp = (a / b).unit();
	Pos3D X = a.unit();//X-axis
	Pos3D Y = (perp / a).unit();//Y-axis
	return angle(X, Y, b);
}
bool inner_check(const Pos3D& d1, const Pos3D& d2, const Pos3D& t, const Pos3D& nrm) {
	return ccw(O3D, d1, t, nrm) >= 0 && ccw(O3D, d2, t, nrm) <= 0;
}
bool connectable(const Polyhedron& P, const Pos3D& a, const Pos3D& b, const int& i, const int& j, const ld& R = 0) {
	if (a == b) return 1;
	Pos3D perp = (a / b).unit();
	Polyhedron inxs;
	for (int k = 0; k < N; k++) {
		if (k == i || k == j) continue;
		Pos3D axis = (P[k] / perp).unit();
		if (zero(axis.Euc())) {
			if ((ld)P[k].r / R < PI * .5) continue;
			else return 0;
		}
		bool f = plane_circle_intersection(perp, P[k], inxs);
		if (!f) continue;
		if (inxs.size() == 1) {
			if (inner_check(a, b, inxs[0], perp)) return 0;
			continue;
		}
		Pos3D hi = inxs[0], lo = inxs[1];
		Pos3D mid = (perp / axis).unit();
		if (ccw(O3D, mid, lo, perp) < 0) std::swap(hi, lo);
		if (inner_check(lo, hi, a, perp) || inner_check(lo, hi, b, perp)) return 0;
		if (inner_check(a, b, lo, perp) || inner_check(a, b, hi, perp)) return 0;
	}
	return 1;
}
Polyhedron circle_circle_intersections(const Pos3D& p, const Pos3D& q, const ld& R = 0) {
	ld ang = angle(p, q);
	ld t1 = (ld)p.r / R;
	ld t2 = (ld)q.r / R;
	assert(t1 <= PI * .5); assert(t2 <= PI * .5);
	if (t1 + t2 < ang || std::abs(t1 - t2) >= ang) return {};
	ld d1 = cos(t1);
	ld d2 = cos(t2);
	Plane p1 = plane(p, -d1);
	Plane p2 = plane(q, -d2);
	Line3D inx;
	if (intersection(p1, p2, inx) != 1) return {};
	ld w = inx.p0.mag();
	ld h = sqrt(1 - sq(w));
	Pos3D v = inx.dir;
	return { inx.p0 + v * h, inx.p0 - v * h };
}
ld spherical_pythagorean(const ld& a, const ld& b, const ld& A, const ld& B) {
	assert(sin(a) > 0); assert(cos(b) > 0);
	ld cosc = fit(cos(a) / cos(b), -1, 1);
	ld c = acos(cosc);
	ld sinC = sin(c) * sin(A) / sin(a);
	return asin(sinC);
}
bool inner_check(const Pos3D& p, const Pos3D& q, const ld& R = 0) {
	ld ang = angle(p, q);
	ld t1 = (ld)p.r / R;
	ld t2 = (ld)q.r / R;
	assert(t1 <= PI * .5); assert(t2 <= PI * .5);
	if (p.r > q.r && std::abs(t1 - t2) >= ang) return 1;
	return 0;
}
Polyhedron tangents(const Pos3D& p, const Pos3D& q, const ld& R = 0) {
	if (zero((p / q).Euc())) return {};
	Polyhedron ret;
	if (!p.r) {
		ld a = angle(q, p);
		ld b = (ld)q.r / R;
		if (a + b >= PI) return {};
		if (a > PI * .5) {
			Pos3D pp = -p;
			Polyhedron tmp = tangents(pp, q);
			return { tmp[1], tmp[0] };
		}
		ld A = PI * .5;
		ld B = asin(fit(sin(A) * sin(b) / sin(a), -1, 1));
		ld C = std::abs(spherical_pythagorean(a, b, A, B));
		Pos3D perp = (q / p).unit();
		Pos3D m = q.rodrigues_rotate(b, perp);
		Pos3D hi = m.rodrigues_rotate(C, q);
		Pos3D lo = m.rodrigues_rotate(-C, q);
		return { hi, lo };
	}
	assert(p.r && q.r);
	ld ang = angle(p, q);
	ld t1 = (ld)p.r / R;
	ld t2 = (ld)q.r / R;
	ld ful = std::abs(ang) + std::abs(t1) + std::abs(t2);
	Polyhedron ptan, qtan, tmp;
	Pos3D rq = -q;
	rq.r = q.r;
	bool f = inner_check(p, rq) || inner_check(rq, p);
	if (ang > (t1 + t2) && !f) {//inner tangent
		if (p.r == q.r) {
			Pos3D m = ((p + q) * .5).unit(); m.r = 0;
			ptan = tangents(m, p);
			qtan = tangents(m, q);
			if (ptan.size() == 2 && qtan.size() == 2) {
				tmp = { ptan[0], qtan[0], ptan[1], qtan[1] };
				ret.insert(ret.end(), tmp.begin(), tmp.end());
			}
		}
		else {
			ld a_ = angle(p, q);
			ld bp = (ld)p.r / R;
			ld bq = (ld)q.r / R;
			ld ap = atan2(sin(a_), (sin(bq) / sin(bp)) + cos(a_));
			Pos3D perp = (p / q).unit();
			Pos3D m = p.rodrigues_rotate(ap, perp); m.r = 0;
			ptan = tangents(m, p);
			qtan = tangents(m, q);
			if (ptan.size() == 2 && qtan.size() == 2) {
				tmp = { ptan[0], qtan[0], ptan[1], qtan[1] };
				ret.insert(ret.end(), tmp.begin(), tmp.end());
			}
		}
	}
	if (ful < PI) {//outer tangent
		if (p.r == q.r) {
			Pos3D top = (q - p).unit(); top.r = 0;
			qtan = tangents(top, q);
			ptan = tangents(top, p);
			if (ptan.size() == 2 && qtan.size() == 2) {
				tmp = { ptan[0], qtan[0], ptan[1], qtan[1] };
				ret.insert(ret.end(), tmp.begin(), tmp.end());
			}
		}
		else {
			Pos3D L = p.r > q.r ? p : q;
			Pos3D S = p.r > q.r ? q : p;
			ld a_ = angle(L, S);
			ld bl = (ld)L.r / R;
			ld bs = (ld)S.r / R;
			ld as = atan2(sin(a_), (sin(bl) / sin(bs)) - cos(a_));
			Pos3D perp = (L / S).unit();
			Pos3D top = S.rodrigues_rotate(as, perp);
			top.r = 0;
			Polyhedron ltan = tangents(top, L);
			Polyhedron stan = tangents(top, S);
			if (ltan.size() == 2 && stan.size() == 2) {
				if (p.r > q.r) tmp = { ltan[0], stan[0], ltan[1], stan[1] };
				else if (p.r < q.r) tmp = { stan[0], ltan[0], stan[1], ltan[1] };
				ret.insert(ret.end(), tmp.begin(), tmp.end());
			}
		}
	}
	return ret;
}

struct Planar {
	Pos3D norm, p0;
	Planar(Pos3D n = Pos3D(0, 0, 0), Pos3D p0_ = Pos3D(0, 0, 0)) : norm(n), p0(p0_) {}
};
Planar P(const Pos3D& p1, const Pos3D& p2, const Pos3D& p3) {
	Pos3D norm = (p2 - p1) / (p3 - p2);
	return Planar(norm, p1);
}
Planar P(std::vector<Pos3D>& tri) {
	Pos3D p1 = tri[0], p2 = tri[1], p3 = tri[2];
	Pos3D norm = (p2 - p1) / (p3 - p2);
	return Planar(norm, p1);
}
Pos3D intersection(const Planar& S, const Line3D& l) {
	ld det = S.norm * l.dir;
	if (zero(det)) return { INF, INF, INF };
	ld t = (S.norm * S.p0 - S.norm * l.p0) / det;
	return l.p0 + (l.dir * t);
}
struct Face {//refer to BIGINTEGER
	int v[3];
	Face(int a, int b, int c) { v[0] = a; v[1] = b; v[2] = c; }
	Pos3D norm(const std::vector<Pos3D>& P) const {
		return cross(P[v[0]], P[v[1]], P[v[2]]);
	}
	bool visible(const std::vector<Pos3D>& P, int i) const {
		return (P[i] - P[v[0]]) * norm(P) > 0;
	}
};
std::vector<Face> H3D;
std::vector<Face> ConvexHull3D(const std::vector<Pos3D>& P) {//refer to BIGINTEGER
	int sz = P.size();
	std::vector<std::vector<int>> vis(sz);
	for (int i = 0; i < sz; i++) vis[i].resize(sz);
	std::vector<Face> cur;
	cur.push_back(Face(0, 1, 2));
	cur.push_back(Face(2, 1, 0));
	for (int i = 3; i < sz; i++) {
		std::vector<Face> next;
		for (int j = 0; j < cur.size(); j++) {
			Face& f = cur[j];
			int ret = f.visible(P, i);
			if (!ret) next.push_back(f);
			for (int k = 0; k < 3; k++) vis[f.v[k]][f.v[(k + 1) % 3]] = ret;
		}
		for (int j = 0; j < cur.size(); j++) {
			for (int k = 0; k < 3; k++) {
				int a = cur[j].v[k], b = cur[j].v[(k + 1) % 3];
				if (vis[a][b] != vis[b][a] && vis[a][b])
					next.push_back(Face(a, b, i));
			}
		}
		cur = next;
	}
	return cur;
}
ld sc[4];
void get_angle(ld sc[], const Pos3D& norm) {
	ld a = norm.x, b = norm.y, c = norm.z;
	ld angle1 = -atan2(b, a);
	ld dx = sqrtl(a * a + b * b);
	ld angle2 = -atan2(dx, c);
	sc[0] = sin(angle1);
	sc[1] = cos(angle1);
	sc[2] = sin(angle2);
	sc[3] = cos(angle2);
	return;
}
void update_sc(const Pos3D& p) {
	ld angle1 = -atan2l(p.y, p.x);
	ld dx = sqrtl(p.x * p.x + p.y * p.y);
	ld angle2 = -atan2l(dx, p.z);
	sc[0] = sinl(angle1);
	sc[1] = cosl(angle1);
	sc[2] = sinl(angle2);
	sc[3] = cosl(angle2);
	return;
}
Pos3D rotate(ld sc[], const Pos3D& p) {//project to xy_plane
	ld x = p.x * sc[1] - p.y * sc[0], y = p.x * sc[0] + p.y * sc[1], z = p.z;
	return Pos3D(z * sc[2] + x * sc[3], y, z * sc[3] - x * sc[2]);
}
void rotate(ld sc[], const Pos3D& p, std::vector<Pos>& C) {//project to xy_plane
	ld x = p.x * sc[1] - p.y * sc[0], y = p.x * sc[0] + p.y * sc[1], z = p.z;
	Pos3D q = Pos3D(z * sc[2] + x * sc[3], y, z * sc[3] - x * sc[2]);
	C.push_back(Pos(q.x, q.y));
	return;
}
Pos3D rotate(const Pos3D& p) {
	ld x = p.x * sc[1] - p.y * sc[0], y = p.x * sc[0] + p.y * sc[1], z = p.z;
	return Pos3D(z * sc[2] + x * sc[3], y, z * sc[3] - x * sc[2]);
}
Pos projecting2D(ld sc[], const Pos3D& p) {//project to xy_plane
	ld x = p.x * sc[1] - p.y * sc[0], y = p.x * sc[0] + p.y * sc[1], z = p.z;
	Pos3D q = Pos3D(z * sc[2] + x * sc[3], y, z * sc[3] - x * sc[2]);
	return Pos(q.x, q.y);
}
Pos convert(Pos3D p, const Pos3D& v) { p = rotate(p - v); return Pos(p.x, p.y); }
Pos convert(Pos3D q, const Pos3D& p, const Pos3D& v) { update_sc(p); return convert(q, v); }
Pos3D recover(const Pos& p2D, const Pos3D& v) {
	ld x = p2D.x * -sc[3];
	ld y = p2D.y;
	ld z = p2D.x * sc[2];
	Pos3D p = Pos3D(x * -sc[1] + y * sc[0], x * sc[0] + y * sc[1], z);
	return p + v;
}