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
#include <set>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
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

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

//#define DEBUG
#ifdef DEBUG
#endif

int N, M, K, Q;
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
const Polygon T[7] = {
	{ Pos(0, 0), Pos(100, 0), Pos(50, 50) },
	{ Pos(0, 0), Pos(100, 0), Pos(50, 50) },
	{ Pos(0, 0), Pos(50, 0), Pos(0, 50) },
	{ Pos(25, 0), Pos(50, 25), Pos(25, 50), Pos(0, 25) },
	{ Pos(0, 0), Pos(50, 0), Pos(75, 25), Pos(25, 25) },
	{ Pos(0, 0), Pos(50, 0), Pos(25, 25) },
	{ Pos(0, 0), Pos(50, 0), Pos(25, 25) }
};
// TR0 | TR0 | TR1 | SQ | TZ | TR2 | TR2
#define SQ 3
#define TZ 4
ld a_ = 625;
const ld A[7] = { a_ * 4, a_ * 4, a_ * 2, a_ * 2, a_ * 2, a_, a_ };
const Polygon TZ2 = { Pos(25, 0), Pos(75, 0), Pos(50, 25), Pos(0, 25) };
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool operator == (const Polygon& P, const Polygon& Q) {
	if (P.size() != Q.size()) return 0;
	int sz = P.size();
	for (int i = 0; i < sz; i++) if (P[i] != Q[i]) return 0;
	return 1;
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
//ld dist(const Pos& d1, const Pos& d2, const Pos& t, bool f = 0) {
//	if (!f) return cross(d1, d2, t) / (d1 - d2).mag();
//	if (sign(projection(d1, d2, d2, t)) <= 0 &&
//		sign(projection(d2, d1, d1, t)) <= 0)
//		return std::abs(cross(d1, d2, t)) / (d1 - d2).mag();
//	return std::min((d1 - t).mag(), (d2 - t).mag());
//}
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);return (p1 * a2 + p2 * a1) / (a1 + a2); }
//ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
ld rad(const Pos& p1, const Pos& p2, const Pos& p3) { return rad(p3 - p2, p1 - p2); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
ld area(const Polygon& H) {
	int sz = H.size();
	ld a = 0;
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a * .5;
}
//int inner_check(Pos H[], const int& sz, const Pos& p) {//concave
//	int cnt = 0;
//	for (int i = 0; i < sz; i++) {
//		Pos cur = H[i], nxt = H[(i + 1) % sz];
//		if (on_seg_strong(cur, nxt, p)) return 1;
//		if (zero(cur.y - nxt.y)) continue;
//		if (nxt.y < cur.y) std::swap(cur, nxt);
//		if (nxt.y - TOL < p.y || cur.y > p.y) continue;
//		cnt += ccw(cur, nxt, p) > 0;
//	}
//	return (cnt & 1) * 2;
//}
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
//Polygon graham_scan(Polygon& C) {
//	Polygon H;
//	if (C.size() < 3) {
//		std::sort(C.begin(), C.end());
//		return C;
//	}
//	std::swap(C[0], *min_element(C.begin(), C.end()));
//	std::sort(C.begin() + 1, C.end(), [&](const Pos& p, const Pos& q) -> bool {
//		int ret = ccw(C[0], p, q);
//		if (!ret) return (C[0] - p).Euc() < (C[0] - q).Euc();
//		return ret > 0;
//		}
//	);
//	C.erase(unique(C.begin(), C.end()), C.end());
//	int sz = C.size();
//	for (int i = 0; i < sz; i++) {
//		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i]) <= 0)
//			H.pop_back();
//		H.push_back(C[i]);
//	}
//	return H;
//}
Polygon convex_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
	Polygon qs;
	int n = ps.size();
	for (int i = 0; i < n; i++) {
		Pos p1 = ps[i], p2 = ps[(i + 1) % n];
		int d1 = ccw(b1, b2, p1), d2 = ccw(b1, b2, p2);
		if (d1 >= 0) qs.push_back(p1);
		if (d1* d2 < 0) qs.push_back(intersection(p1, p2, b1, b2));
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
	//Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
	Pos p(const ld& rt = .5) const { return s + dir * rt; }
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
//Segs half_plane_intersection(Segs& HP, const bool& srt = 1) {
//	auto check = [&](Seg& u, Seg& v, Seg& w) -> bool {
//		return w.inner(intersection_(u, v));
//		};
//	if (srt) std::sort(HP.begin(), HP.end());
//	std::deque<Seg> dq;
//	int sz = HP.size();
//	for (int i = 0; i < sz; ++i) {
//		if (i && same_dir(HP[i], HP[(i - 1) % sz])) continue;
//		while (dq.size() > 1 && !check(dq[dq.size() - 2], dq[dq.size() - 1], HP[i])) dq.pop_back();
//		while (dq.size() > 1 && !check(dq[1], dq[0], HP[i])) dq.pop_front();
//		dq.push_back(HP[i]);
//	}
//	while (dq.size() > 2 && !check(dq[dq.size() - 2], dq[dq.size() - 1], dq[0])) dq.pop_back();
//	while (dq.size() > 2 && !check(dq[1], dq[0], dq[dq.size() - 1])) dq.pop_front();
//	sz = dq.size();
//	if (sz < 3) return {};
//	std::vector<Seg> HPI;
//	for (int i = 0; i < sz; ++i) HPI.push_back(dq[i]);
//	return HPI;
//}
//Segs half_plane_intersection(const Segs& P, const Segs& Q) {
//	Segs HP(P.size() + Q.size());
//	std::merge(P.begin(), P.end(), Q.begin(), Q.end(), HP.begin());
//	return half_plane_intersection(HP, 0);
//}
bool seg_overlap(const Seg& p, const Seg& q) {
	bool f0 = collinear(p.s, p.e, q.s, q.e);
	bool f1 = sign(dot(p.s, p.e, q.s, q.e)) > 0;
	//bool f2 = on_seg_strong(p.s, p.e, q.s) || on_seg_strong(p.s, p.e, q.e);
	return f0 && f1;
}
bool cut_seg(const Seg& se, const Polygon& C, Segs& ret) {
	int szc = C.size();
	ld s = 0, e = 1;
	Polygon V;
	for (int i = 0; i < szc; i++) {
		int i0 = i, i1 = (i + 1) % szc;
		Seg q = Seg(C[i0], C[i1]);
		if (seg_overlap(se, q)) {
			s = fit(projection(se.s, se.e, se.s, C[i0]), 0, 1);
			e = fit(projection(se.s, se.e, se.s, C[i1]), 0, 1);
			V.push_back(Pos(s, e));
		}
		//		if (ccw(se.s, se.e, C[i0], C[i1])) {
		////#ifdef DEBUG
		////			std::cout << "se.s:: " << se.s << " ";
		////			std::cout << "se.e:: " << se.e << " ";
		////			std::cout << "C[i0]:: " << C[i0] << " ";
		////			std::cout << "C[i1]:: " << C[i1] << "\n";
		////#endif
		//			if (on_seg_weak(se.s, se.e, C[i0]) ||
		//				on_seg_weak(se.s, se.e, C[i1])) {
		//				std::cout << "se.s:: " << se.s << "\n";
		//				std::cout << "se.e:: " << se.e << "\n";
		//				std::cout << "C[i0]:: " << C[i0] << "\n";
		//				std::cout << "C[i1]:: " << C[i1] << "\n";
		//				return 0;
		//			}
		//			if (on_seg_weak(C[i0], C[i1], se.s) ||
		//				on_seg_weak(C[i0], C[i1], se.e)) {
		//				std::cout << "se.s:: " << se.s << "\n";
		//				std::cout << "se.e:: " << se.e << "\n";
		//				std::cout << "C[i0]:: " << C[i0] << "\n";
		//				std::cout << "C[i1]:: " << C[i1] << "\n";
		//				return 0;
		//			}
		//		}
	}
	std::sort(V.begin(), V.end());
	V.push_back(Pos(1, 1));
	ld hi = 0;
	for (const Pos& p : V) {
		if (hi < p.LO) ret.push_back(Seg(se.p(hi), se.p(p.LO))), hi = p.HI;
		else hi = std::max(hi, p.HI);
	}
	return 1;
}
bool substract(const Polygon& P, const Polygon& C, Polygon& ret) {
	int szp = P.size();
	int szc = C.size();
	Segs VP, VC;
#ifdef DEBUG
	std::cout << "sub start::\n";
#endif
	for (int i = 0; i < szp; i++) {
		int i0 = i, i1 = (i + 1) % szp;
		Seg se = Seg(P[i0], P[i1]);
		Segs vse;
		bool f = cut_seg(se, C, vse);
#ifdef DEBUG
		std::cout << "f0:: " << f << "\n";
#endif
		if (!f) return 0;
		for (const Seg& v : vse) {
			if (v.s == v.e) continue;
			VP.push_back(v);
		}
	}
	for (int i = 0; i < szc; i++) {
		int i0 = i, i1 = (i + 1) % szc;
		Seg se = Seg(C[i0], C[i1]);
		Segs vse;
		bool f = cut_seg(se, P, vse);
#ifdef DEBUG
		std::cout << "f1:: " << f << "\n";
#endif
		if (!f) return 0;
		for (const Seg& v : vse) {
			if (v.s == v.e) continue;
			VC.push_back(v);
		}
	}
#ifdef DEBUG
	std::cout << "seg done::\n";
#endif
	int p = VP.size(), c = VC.size();
	for (int i = 0; i < p; i++) {
		if (VP[(i - 1 + p) % p].e != VP[i].s) {
			std::rotate(VP.begin(), VP.begin() + i, VP.end());
			break;
		}
	}
#ifdef DEBUG
	std::cout << "VP::\n";
	for (const Seg& se : VP) {
		std::cout << "  [(" << se.s.x << ", " << se.s.y << "), ";
		std::cout << "(" << se.e.x << ", " << se.e.y << ")]\n";
	}
#endif
	for (int i = 0; i < p - 1; i++) if (VP[i].e != VP[i + 1].s) return 0;
#ifdef DEBUG
	std::cout << "P done::\n";
#endif
	std::reverse(VC.begin(), VC.end());
	for (Seg& se : VC) std::swap(se.s, se.e);
	for (int i = 0; i < c; i++) {
		if (VC[(i - 1 + c) % c].e != VC[i].s) {
			std::rotate(VC.begin(), VC.begin() + i, VC.end());
			break;
		}
	}
#ifdef DEBUG
	std::cout << "VC::\n";
	for (const Seg& se : VC) {
		std::cout << "  [(" << se.s.x << ", " << se.s.y << "), ";
		std::cout << "(" << se.e.x << ", " << se.e.y << ")]\n";
	}
#endif
	for (int i = 0; i < c - 1; i++) if (VC[i].e != VC[i + 1].s) return 0;
#ifdef DEBUG
	std::cout << "C done::\n";
#endif
	assert(VP[0].s == VC.back().e);
	assert(VC[0].s == VP.back().e);
	for (const Seg& se : VP) ret.push_back(se.s);
	for (const Seg& se : VC) ret.push_back(se.s);
	Polygon tmp = ret;
	std::sort(ret.begin(), ret.end());
	int sz = ret.size();
	for (int i = 0; i < sz - 1; i++) if (ret[i] == ret[i + 1]) return 0;
#ifdef DEBUG
	std::cout << "dup done::\n";
#endif
	Vbool F(sz, 0);
	for (int i = 0; i < sz; i++) {
		int i0 = (i - 1 + sz) % sz, i1 = i, i2 = (i + 1) % sz;
		if (!ccw(tmp[i0], tmp[i1], tmp[i2])) F[i1] = 1;
	}
	ret.clear();
	for (int i = 0; i < sz; i++) if (!F[i]) ret.push_back(tmp[i]);
#ifdef DEBUG
	std::cout << "prune done::\n";
#endif
	return 1;
}
Polygon rotate_and_move(const Polygon& P, const int& i, const ld& t, const Pos& pv) {
	Polygon C = P;
	//Pos vec = P[i];
	for (Pos& p : C) p = p.rot(t);
	Pos vec = C[i] - pv;
	for (Pos& p : C) p -= vec;
	return C;
}
bool two_polygon_cmp(const Polygon& P, const int& s) {
	int sz = P.size();
	const Polygon& Q = T[s];
	if (Q.size() != sz) return 0;
	if (s == SQ) {//sqaure check
		assert(sz == 4);
		ld lq = (Q[0] - Q[1]).mag();
		for (int i = 0; i < 4; i++) {
			if (sign(dot(P[i], P[(i + 1) % sz], P[(i + 2) % sz]))) return 0;
			ld lp = (P[i] - P[(i + 1) % sz]).mag();
			if (!eq(lp, lq)) return 0;
		}
		return 1;
	}
	if (s == TZ) {//parallelogram check
		assert(sz == 4);
		const Pos& p0 = P[0], & p1 = P[1];
		const Pos& p2 = P[2], & p3 = P[3];
		ld d0 = (Q[0] - Q[1]).mag();
		ld d1 = (Q[1] - Q[2]).mag();
		if (d1 < d0) std::swap(d0, d1);
		if (ccw(p0, p1, p3, p2)) return 0;
		if (ccw(p1, p2, p0, p3)) return 0;
		ld t = norm(rad(p0, p1, p2));
		if (!eq(t, PI * .25) && !eq(t, PI * .75)) return 0;
		ld l0 = (p0 - p1).mag();
		ld l1 = (p1 - p2).mag();
		if (l1 < l0) std::swap(l0, l1);
		if (eq(l0, d0) && eq(l1, d1)) return 1;
		return 0;
	}
	//triangle checkk
	assert(sz == 3);
	Vld vp, vq;
	for (int i = 0; i < 3; i++) {
		ld lp = (P[i] - P[(i + 1) % sz]).mag();
		vp.push_back(lp);
		ld lq = (Q[i] - Q[(i + 1) % sz]).mag();
		vq.push_back(lq);
	}
	std::sort(vp.begin(), vp.end());
	std::sort(vq.begin(), vq.end());
	for (int i = 0; i < 3; i++) if (!eq(vp[i], vq[i])) return 0;
	return 1;
}
std::vector<Polygon> candidates(const int& x, const Polygon& P, const int& i) {
	const Polygon& C = T[x];
	std::vector<Polygon> RET;
	int sz = P.size();
	int i0 = (i - 1 + sz) % sz, i1 = i, i2 = (i + 1) % sz;
	int tq = ccw(P[i0], P[i1], P[i2]);
	int fc = sign(dot(P[i0], P[i1], P[i2]));
	if (tq <= 0) return {};
	Pos u, v;
	ld t, a;
	Polygon R;
#ifdef DEBUG
	//std::cout << "candidates::\n";
	//std::cout << "x:: " << x << "\n";
	//std::cout << "i:: " << i << "\n";
	//std::cout << "P::\n";
	//for (const Pos& p : P) {
	//	std::cout << "  (" << p.x << ", " << p.y << "),\n";
	//}
	//std::cout << "C::\n";
	//for (const Pos& p : C) {
	//	std::cout << "  (" << p.x << ", " << p.y << "),\n";
	//}
#endif
	if (x == SQ) {
		if (fc < 0) return {};
		if (!fc) {
			u = C[1] - C[0];
			v = P[i2] - P[i1];
			t = rad(u, v);
			R = rotate_and_move(C, 0, t, P[i1]);
			a = area(sutherland_hodgman(P, R));
			if (eq(a, A[x])) RET.push_back(R);
		}
		else {
			u = C[1] - C[0];
			v = P[i2] - P[i1];
			t = rad(u, v);
			R = rotate_and_move(C, 0, t, P[i1]);
			a = area(sutherland_hodgman(P, R));
			if (eq(a, A[x])) RET.push_back(R);
			u = C[3] - C[0];
			v = P[i0] - P[i1];
			t = rad(u, v);
			R = rotate_and_move(C, 0, t, P[i1]);
			a = area(sutherland_hodgman(P, R));
			if (eq(a, A[x])) RET.push_back(R);
		}
	}
	else if (x == TZ) {
		for (int j = 0; j < 2; j++) {
			for (const Polygon& Z : { C, TZ2 }) {
				u = Z[j + 1] - Z[j];
				v = P[i0] - P[i1];
				t = rad(u, v);
				R = rotate_and_move(Z, j, t, P[i1]);
				a = area(sutherland_hodgman(P, R));
				if (eq(a, A[x])) RET.push_back(R);
				u = Z[(j + 3) % 4] - Z[j];
				v = P[i2] - P[i1];
				t = rad(u, v);
				R = rotate_and_move(Z, j, t, P[i1]);
				a = area(sutherland_hodgman(P, R));
				if (eq(a, A[x])) RET.push_back(R);
			}
		}
	}
	else {//TR
		//#ifdef DEBUG
		//		std::cout << "triangle::\n";
		//#endif
		for (int j = 0; j < 3; j++) {
			u = C[(j + 1) % 3] - C[j];
			v = P[i2] - P[i1];
			t = rad(u, v);
			R = rotate_and_move(C, j, t, P[i1]);
			//#ifdef DEBUG
			//			std::cout << "t:: " << t << "\n";
			//			std::cout << "C[j]:: " << C[j] << "\n";
			//			std::cout << "C[j + 1]:: " << C[(j + 1) % 3] << "\n";
			//			std::cout << "P[i1]:: " << P[i1] << "\n";
			//			std::cout << "P[i2]:: " << P[i2] << "\n";
			//			std::cout << "R::\n";
			//			for (const Pos& p : R) {
			//				std::cout << "  (" << p.x << ", " << p.y << "),\n";
			//			}
			//#endif
			a = area(sutherland_hodgman(P, R));
			if (eq(a, A[x])) RET.push_back(R);
			u = C[(j + 2) % 3] - C[j];
			v = P[i0] - P[i1];
			t = rad(u, v);
			R = rotate_and_move(C, j, t, P[i1]);
			//#ifdef DEBUG
			//			std::cout << "t:: " << t << "\n";
			//			std::cout << "C[j]:: " << C[j] << "\n";
			//			std::cout << "C[j + 2]:: " << C[(j + 2) % 3] << "\n";
			//			std::cout << "P[i1]:: " << P[i1] << "\n";
			//			std::cout << "P[i0]:: " << P[i0] << "\n";
			//			std::cout << "R::\n";
			//			for (const Pos& p : R) {
			//				std::cout << "  (" << p.x << ", " << p.y << "),\n";
			//			}
			//#endif
			a = area(sutherland_hodgman(P, R));
			if (eq(a, A[x])) RET.push_back(R);
		}
		//#ifdef DEBUG
		//		std::cout << "triangle done::\n";
		//		for (const Polygon& C : RET) {
		//			std::cout << "C::\n";
		//			for (const Pos& p : C) {
		//				std::cout << "  (" << p.x << ", " << p.y << "),\n";
		//			}
		//		}
		//#endif
	}
	sz = RET.size();
	//#ifdef DEBUG
	//	std::cout << "sz::" << sz << "\n";
	//#endif
	Vbool F(sz, 0);
	for (int j = 0; j < sz; j++) {
		for (int k = j + 1; k < sz; k++) {
			if (RET[j] == RET[k]) {
				F[k] = 1;
			}
		}
	}
	//#ifdef DEBUG
	//	std::cout << "compare done\n";
	//#endif
	std::vector<Polygon> TMP;
	for (int j = 0; j < sz; j++) if (!F[j]) TMP.push_back(RET[j]);
	return TMP;
}
std::set<int> I;
#define POLYGON_DEBUG
#ifdef POLYGON_DEBUG
Polygon S[7];
#endif
bool dfs(int d, const Polygon& P) {
	if (d == 6) {
		assert(I.size() == 6);
		for (int i = 0; i < 7; i++) {
			if (I.count(i)) continue;
			if (two_polygon_cmp(P, i)) {
#ifdef DEBUG
				std::cout << "d:: " << d << "\n";
				std::cout << "final floor\n";
#endif
#ifdef POLYGON_DEBUG
				S[d] = P;
				for (int s = 0; s < 7; s++) {
					std::cout << "S" << s << " = [\n";
					for (Pos& p : S[s]) {
						std::cout << "  (" << p.x << ", " << p.y << "),\n";
					}
					std::cout << "]\n";
				}
#endif
				return 1;
			}
			break;
		}
		return 0;
	}
	//std::cout << "\nd:: " << d << "\n";
#ifdef DEBUG
#endif
	for (int i = 0; i < 7; i++) {
		if (I.count(i)) continue;
		//std::cout << "i:: " << i << "\n";
		I.insert(i);
		int sz = P.size();
#ifdef DEBUG
		std::cout << "sz:: " << sz << "\n";
#endif
		for (int j = 0; j < sz; j++) {
			std::vector<Polygon> CC = candidates(i, P, j);
#ifdef DEBUG
			std::cout << "CC.sz:: " << CC.size() << "\n";
#endif
			for (const Polygon& C : CC) {
#ifdef DEBUG
				std::cout << "C:: \n";
				std::cout << "C.sz:: " << C.size() << "\n";
				for (const Pos& p : C) {
					std::cout << "  (" << p.x << ", " << p.y << "),\n";
				}
#endif
				Polygon R;
				bool f = substract(P, C, R);
#ifdef DEBUG
				std::cout << "f:: " << f << "\n";
#endif
				//#ifdef POLYGON_DEBUG
				//				S[d] = C;
				//#endif
				//#ifdef POLYGON_DEBUG
				//				if (!f) {
				//					S[d] = C;
				//					for (int s = 0; s <= d; s++) {
				//						std::cout << "S" << s << " = [\n";
				//						for (Pos& p : S[s]) {
				//							std::cout << "  (" << p.x << ", " << p.y << "),\n";
				//						}
				//						std::cout << "]\n";
				//					}
				//					std::cout << "P = [\n";
				//					for (const Pos& p : P) {
				//						std::cout << "  (" << p.x << ", " << p.y << "),\n";
				//					}
				//					std::cout << "]\n\n";
				//				}
				//#endif
				if (!f) continue;
#ifdef POLYGON_DEBUG
				S[d] = C;
				std::cout << "S" << d << " = [\n";
				for (Pos& p : S[d]) {
					std::cout << "  (" << p.x << ", " << p.y << "),\n";
				}
				std::cout << "]\n";
				std::cout << "R = [\n";
				for (Pos& p : R) {
					std::cout << "  (" << p.x << ", " << p.y << "),\n";
				}
				std::cout << "]\n\n";
#endif
#ifdef DEBUG
				std::cout << "d:: " << d << " go next floor\n";
				for (Pos& p : S[d]) {
					std::cout << "  (" << p.x << ", " << p.y << "),\n";
				}
#endif
				f = dfs(d + 1, R);
				if (f) return 1;
			}
		}
		I.erase(i);
	}
	return 0;
}
bool tangram_solver(const Polygon& P) {
	ld A = area(P); if (!eq(A, 1e4)) return 0;
	int sz = P.size(); if (sz > 23) return 0;
	return dfs(0, P);
}//badk tracking
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	//freopen("..//", "r", stdin);
	//freopen("../../tests/tangram_out.txt", "w", stdout);
	std::cin >> N;
	Polygon P(N);
	for (Pos& p : P) {
		int a, b, c, d;
		std::cin >> a >> b >> c >> d;
		p.x = a + b * sqrt(2);
		p.y = c + d * sqrt(2);
	}
	ld A = area(P);
	if (A < 0) std::reverse(P.begin(), P.end());
	bool f = tangram_solver(P);
	std::cout << (f ? "YES\n" : "NO\n");
	return;
}
int main() { solve(); return 0; }//boj22654 Tangram