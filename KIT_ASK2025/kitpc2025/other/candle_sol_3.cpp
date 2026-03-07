#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <iomanip>
#include <chrono>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<size_t> Vidx;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
const ld INF = 1e17;
const ld TOL = 1e-15;
const int LEN = 150;
const ld PI = acosl(-1);
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ll sq(const ll& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo = 0, const ld& hi = 1) { return std::min(hi, std::max(lo, x)); }

#define STRONG 0
#define WEAK 1

#define LO x
#define HI y

#define BLACK   0
#define RED     (1 << 0)
#define GREEN   (1 << 1)
#define BLUE    (1 << 2)
#define YELLOW  (RED | GREEN)
#define MAGENTA (RED | BLUE)
#define CYAN    (GREEN | BLUE)
#define WHITE   (RED | GREEN | BLUE)

//=======================DEBUG MACRO=======================//
//#define DEBUG
#ifdef DEBUG
#endif

//#define ASSERT
#ifdef ASSERT
#endif

//#define POLYGON_CHECK

#define NAIVE
#ifndef NAIVE
#define FAST
#endif
//=======================DEBUG MACRO=======================//

int N, M, K, Q;
ld A[1 << 3];
int I[1 << 3];
int C[LEN];
struct Pos {
	ld x, y;
	Pos(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& n) const { return { x * n, y * n }; }
	Pos operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pos& p) const { return { x * p.x + y * p.y }; }
	ld operator / (const Pos& p) const { return { x * p.y - y * p.x }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { -x, -y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& scale) { x *= scale; y *= scale; return *this; }
	Pos& operator /= (const ld& scale) { x /= scale; y /= scale; return *this; }
	Pos rot(const ld& t) { return { x * cosl(t) - y * sinl(t), x * sinl(t) + y * cosl(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrtl(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2l(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "  (" << p.x << ", " << p.y << ")"; return os; }
} L[1 << 3]; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon P[LEN];
Polygon H[LEN];
std::vector<Polygon> T[1 << 3];
Vld R[1 << 3];
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return dot(d1, d2, d1, d3) / (d2 - d1).Euc(); }
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
ld area(const Polygon& H) {
	int sz = H.size();
	ld a = 0;
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a * .5;
}
void norm(Polygon& H) {
	ld a = area(H);
	assert(!zero(a));
	if (a < 0) std::reverse(H.begin(), H.end());
	return;
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
Polygon polygon_cut(const Polygon& ps, const Pos& b1, const Pos& b2) {
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
	Polygon ret = C;
	for (int i = 0; i < sz; i++) {
		Pos b1 = clip[i], b2 = clip[(i + 1) % sz];
		ret = polygon_cut(ret, b1, b2);
	}
	return ret;
}
struct Seg {
	Pos s, e;
	int i;
	Seg(Pos s_ = Pos(), Pos e_ = Pos(), int i_ = -1) : s(s_), e(e_), i(i_) {}
	Pos dir() const { return (s - e).unit(); }
	//bool operator < (const Seg& l) const {
	//	Pos v0 = dir();
	//	Pos v1 = l.dir();
	//	bool f0 = O < v0;
	//	bool f1 = O < v1;
	//	if (f0 != f1) return f0;
	//	if (collinear(s, e, l.s, l.e)) return s == l.s ? e < l.e : s < l.s;
	//	int tq = v0 / v1;
	//	return !tq ? ccw(s, e, l.s) > 0 : tq > 0;
	//}
	bool operator < (const Seg& l) const { return s == l.s ? e < l.e : s < l.s; }
	bool operator == (const Seg& l) const { return s == l.s && e == l.e; }
	Pos p(const ld& rt = .5) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
Seg S[1 << 3][LEN];
typedef std::vector<Seg> Segs;
ld dot(const Seg& p, const Seg& q) { return dot(p.s, p.e, q.s, q.e); }
bool intersect(const Seg& u, const Seg& v) { return intersect(u.s, u.e, v.s, v.e); }
ld intersection(const Seg& s1, const Seg& s2, const bool& f = STRONG) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == WEAK) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
bool inner_check(const Polygon& H, const Pos& q, const Pos& d = Pos(0, 0)) {//convex
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos& p1 = H[i], & p2 = H[(i + 1) % sz];
		if (ccw(p1, p2, q) < 0) return 0;
		if (on_seg_strong(p1, p2, q) && !eq(d.x, d.y)) {
			if (sign((p2 - p1) / d) > 0) return 1;
			else return 0;
		}
	}
	return 1;
}
bool inner_check_concave(const std::vector<Pos>& H, const Pos& p, const Pos& s, const Pos& e) {
	int cnt = 0, sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) {
			//std::cout << "cur:: " << cur << "\n";
			//std::cout << "nxt:: " << nxt << "\n";
			//std::cout << "s:: " << s << "\n";
			//std::cout << "e:: " << e << "\n";
			//assert(collinear(cur, nxt, s, e));
			return dot(cur, nxt, s, e) > 0 ? 1 : 0;
		}
		if (zero(cur.y - nxt.y)) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y - TOL < p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return cnt & 1;
}
Pos get_pos(const Pos& l, const Seg& p, const Seg& q) {
	Pos p1 = p.s, p2 = p.e, q1 = q.s, q2 = q.e;
	if (!inside(p2, l, p1, q1, WEAK) && !inside(p2, l, p1, q2, WEAK)) {
		if (intersect(l, p1, q1, q2) && intersect(l, p2, q1, q2)) return Pos(0, 1);
		else return Pos(0, 0);
	}
	Polygon tri = { p1, p2, l };
	bool in1 = inner_check(tri, q1), in2 = inner_check(tri, q2);
	if (!in1 && !in2) return Pos(0, 0);
	ld r1 = 0, r2 = 1;
	if (in1 && in2) {
		r1 = intersection(p, Seg(l, q1), WEAK);
		r2 = intersection(p, Seg(l, q2), WEAK);
	}
	else if (in1) r1 = intersection(p, Seg(l, q1), WEAK);
	else if (in2) r2 = intersection(p, Seg(l, q2), WEAK);
	else r1 = r2 = 0;
	if (r2 < r1) std::swap(r1, r2);
	return Pos(r1, r2);
}
Vld intersections(const Seg& l, const int& t) {
	int sz = T[t].size();
	Vld ret;
	for (int i = 0; i < sz; i++) {
		const Polygon& Q = T[t][i];
		assert(Q.size() == 3);
		for (int j = 0; j < 3; j++) {
			const Pos& t0 = Q[j], & t1 = Q[(j + 1) % 3];
			//ld x = intersection(l, Seg(t0, t1));
			//if (x > -.5) ret.push_back(x);
			if (collinear(l.s, l.e, t0, t1)) {
				for (const Pos& p : { t0, t1 }) {
					ld ix = projection(l.s, l.e, p);
					ret.push_back(fit(ix));
				}
			}
			else if (ccw(l.s, l.e, t0) * ccw(l.s, l.e, t1) <= 0) {
				ld ix = intersection(l, Seg(t0, t1), WEAK);
				ret.push_back(fit(ix));
			}
		}
	}
	std::sort(ret.begin(), ret.end());
	ret.erase(unique(ret.begin(), ret.end(), eq), ret.end());
	return ret;
}
struct Event {
	ld x;
	int f;
	bool operator < (const Event& o) const { return eq(x, o.x) ? f < o.f : sign(x - o.x) < 0; }
	bool operator == (const Event& o) const { return eq(x, o.x) && f == o.f; }
};
typedef std::vector<Event> Ve;
Pos cen;
Segs Z[1 << 3];
bool cmpt(const Seg& p, const Seg& q) {
	assert(ccw(cen, p.s, p.e));
	assert(ccw(cen, q.s, q.e));
	Pos u = p.s - cen;
	Pos v = q.s - cen;
	return norm(u.rad()) < norm(v.rad());
	//bool f1 = O < u;
	//bool f2 = O < v;
	//if (f1 != f2) return f1;
	//assert(!zero(u / v));
	//return u / v > 0;
}
bool inner_check_bi_search(const int& i, const Seg& se) {
	const Pos& c = L[i], q = se.p();
	ld tm = norm((q - c).rad());
	assert(T[i].size() == R[i].size());
	int sz = T[i].size();
	int s = 0, e = sz - 1, m = 0;
	while (s < e) {
		m = s + e >> 1;
		ld t = R[i][m];
		if (eq(t, tm)) break;
		if (tm < t) e = m;
		else s = m + 1;
	}
	int j = (m - 1 + sz) % sz;
	for (int _ = 0; _ < 3; _++) {
		bool f = inner_check_concave(T[i][j], se.p(), se.s, se.e);
		if (f) return 1;
		j = (j + 1) % sz;
	}
	return 0;
}
void query(const int& q) {
	std::cin >> L[RED] >> L[GREEN] >> L[BLUE];
	memset(A, 0, sizeof A);
	memset(I, -1, sizeof I);
	memset(C, 0, sizeof C);
	for (int c = 1; c < (1 << 3); c++) T[c].clear();
	for (int c = 1; c < (1 << 3); c++) H[c].clear();
	for (int c = 1; c < (1 << 3); c++) Z[c].clear();
	for (int c = 1; c < (1 << 3); c++) R[c].clear();
	for (int c = 1; c < (1 << 3); c <<= 1) {
		int f0 = inner_check(P[0], L[c]), fk = 0;
		if (!f0) continue;
		if (f0) { C[0] |= c; I[c] = 0; }
		for (int k = 1; k <= K; k++) {
			fk = inner_check(P[k], L[c]);
			if (fk) {
				C[0] -= c;
				C[k] |= c;
				I[c] = k;
				break;
			}
		}
		if (fk) continue;
		for (int k = 1; k <= K; k++) {
			const Polygon& H = P[k];
			M = H.size();
			Pos pl = H[0], pr = H[0];
			for (int j = 0; j < M; j++) {
				if (ccw(L[c], pl, H[j]) > 0) pl = H[j];
				if (ccw(L[c], pr, H[j]) < 0) pr = H[j];
			}
			S[c][k] = Seg(pr, pl);
		}
		int sz = P[0].size();
		for (int i = 0; i < sz; i++) {
			const Pos& u = P[0][i], & v = P[0][(i + 1) % sz];
			assert(ccw(L[c], u, v) > 0);
			Seg w = Seg(u, v);
			Polygon VP = { Pos(0, 0) };
			for (int k = 1; k <= K; k++) {
				Pos se = get_pos(L[c], w, S[c][k]);
				if (!eq(se.x, se.y)) VP.push_back(se);
			}
			VP.push_back(Pos(1, 1));
			std::sort(VP.begin(), VP.end());
			ld hi = 0;
			for (const Pos& p : VP) {
				if (hi < p.LO) {
					//if (sign(hi - p.LO) < 0) {
					Pos s = w.p(hi);
					Pos e = w.p(p.LO);
					Z[c].push_back(Seg(s, e));
					hi = p.HI;
				}
				else hi = std::max(hi, p.HI);
			}
		}
		for (int k = 1; k <= K; k++) {
			const Polygon& H = P[k];
			sz = H.size();
			for (int i = 0; i < sz; i++) {
				const Pos& u = H[i], & v = H[(i + 1) % sz];
				if (ccw(L[c], u, v) >= 0) continue;
				Seg w = Seg(v, u);
				Polygon VP = { Pos(0, 0) };
				for (int k_ = 1; k_ <= K; k_++) {
					if (k_ == k) continue;
					Pos se = get_pos(L[c], w, S[c][k_]);
					if (!eq(se.x, se.y)) VP.push_back(se);
				}
				VP.push_back(Pos(1, 1));
				std::sort(VP.begin(), VP.end());
				ld hi = 0;
				for (const Pos& p : VP) {
					if (hi < p.LO) {
						//if (sign(hi - p.LO) < 0) {
						Pos s = w.p(hi);
						Pos e = w.p(p.LO);
						Z[c].push_back(Seg(s, e));
						hi = p.HI;
					}
					else hi = std::max(hi, p.HI);
				}
			}
		}
		cen = L[c]; std::sort(Z[c].begin(), Z[c].end(), cmpt);
		sz = Z[c].size();
		for (int z = 0; z < sz; z++) {
			const Pos& s = Z[c][z].s;
			const Pos& e = Z[c][z].e;
			H[c].push_back(s);
			H[c].push_back(e);
			ld t = norm((s - L[c]).rad());
			T[c].push_back({ L[c], s, e });
			R[c].push_back(t);
		}
		H[c].erase(unique(H[c].begin(), H[c].end()), H[c].end());
		if (H[c][0] == H[c].back()) H[c].pop_back();
		//std::cout << "H["<< c <<"]::\n";
		//for (Pos& p : H[c]) std::cout << p << "\n";
		//std::cout << "H["<< c <<"]::\n";
	}
	if (!I[RED] && !I[GREEN] && !I[BLUE]) {//R & G & B
		Segs VS;
		for (int i = 0; i < 3; i++) {
			int c0 = 1 << i;
			int c1 = (1 << ((i + 1) % 3));
			int c2 = (1 << ((i + 2) % 3));
			int sz = H[c0].size();
			for (int j = 0; j < sz; j++) {
				Vld V = { 0, 1 };
				int j0 = j, j1 = (j + 1) % sz;
				Seg se = Seg(H[c0][j0], H[c0][j1]);
				for (const int& c : { c1, c2 }) {
					Vld tmp = intersections(se, c);
					V.insert(V.end(), tmp.begin(), tmp.end());
				}
				std::sort(V.begin(), V.end());
				V.erase(unique(V.begin(), V.end(), eq), V.end());
				int szv = V.size();
				for (int v = 0; v < szv - 1; v++) {
					ld s = V[v], e = V[v + 1];
					VS.push_back(Seg(se.p(s), se.p(e), c0));
				}
			}
		}
		std::sort(VS.begin(), VS.end());
		VS.erase(unique(VS.begin(), VS.end()), VS.end());
		//std::cout << "VS:: " << VS.size() << "\n";
		for (const Seg& se : VS) {
			bool f = 1;
			for (int i = 0; i < 3; i++) {
				int c = 1 << i;
				if (c == se.i) continue;
				if (!inner_check_bi_search(c, se)) {
					f = 0;
					break;
				}
			}
			//std::cout << "f:: " << f << " ";
			//std::cout << "se:: " << se.s << " " << se.e << "\n";
			if (f) A[WHITE] += se.green();
		}
	}
	for (int i = 0; i < 3; i++) {//R & G, G & B, B & R  3*9*N^2 -> N==200:1,080,000
		int c1 = (1 << ((i + 1) % 3));
		int c2 = (1 << ((i + 2) % 3));
		if (!I[c1] && !I[c2]) {
			int c = c1 | c2;
			Segs VS;
			for (int _ = 0; _ < 2; _++) {
				int sz = H[c1].size();
				for (int j = 0; j < sz; j++) {
					Vld V = { 0, 1 };
					int j0 = j, j1 = (j + 1) % sz;
					Seg se = Seg(H[c1][j0], H[c1][j1]);
					Vld tmp = intersections(se, c2);
					V.insert(V.end(), tmp.begin(), tmp.end());
					std::sort(V.begin(), V.end());
					V.erase(unique(V.begin(), V.end(), eq), V.end());
					int szv = V.size();
					for (int v = 0; v < szv - 1; v++) {
						ld s = V[v], e = V[v + 1];
						VS.push_back(Seg(se.p(s), se.p(e), c1));
					}
				}
				std::swap(c1, c2);
			}
			std::sort(VS.begin(), VS.end());
			VS.erase(unique(VS.begin(), VS.end()), VS.end());
			for (const Seg& se : VS) {
				bool f = 1;
				for (const int& c : { c1, c2 }) {
					if (c == se.i) continue;
					if (!inner_check_bi_search(c, se)) {
						f = 0;
						break;
					}
				}
				if (f) A[c] += se.green();
			}
		}
	}
	for (int c = 1; c < (1 << 3); c <<= 1) {//R, G, B  3*3*N -> N==200:1,800
		if (!I[c]) {
			A[c] += area(H[c]);
			//for (const Polygon& T0 : T[c]) A[c] += area(T0);
		}
	}
	for (int k = 1; k <= K; k++) {
		ld a = area(P[k]);
		int c = C[k];
		A[c] += a;
	}
	if (!I[RED] && !I[GREEN] && !I[BLUE]) {
		for (int c = 1; c < WHITE; c++) A[c] -= A[WHITE];
	}
	for (int i = 0; i < 3; i++) {
		int c1 = (1 << ((i + 1) % 3));
		int c2 = (1 << ((i + 2) % 3));
		if (!I[c1] && !I[c2]) {
			int c = c1 | c2;
			A[c1] -= A[c];
			A[c2] -= A[c];
		}
	}
	A[BLACK] = area(P[0]);
	for (int c = 1; c < (1 << 3); c++) A[BLACK] -= A[c];
	for (int c = 0; c < (1 << 3); c++) A[c] = std::max(A[c], (ld).0);
	//std::cout << "Case #" << q << ":\n";
	//std::cout << "R: " << A[RED] << "\n";
	//std::cout << "G: " << A[GREEN] << "\n";
	//std::cout << "B: " << A[BLUE] << "\n";
	//std::cout << "Y: " << A[YELLOW] << "\n";
	//std::cout << "M: " << A[MAGENTA] << "\n";
	//std::cout << "C: " << A[CYAN] << "\n";
	//std::cout << "W: " << A[WHITE] << "\n";
	//std::cout << "L: " << A[BLACK] << "\n";
	std::cout << A[RED] << "\n";
	std::cout << A[GREEN] << "\n";
	std::cout << A[BLUE] << "\n";
	std::cout << A[YELLOW] << "\n";
	std::cout << A[MAGENTA] << "\n";
	std::cout << A[CYAN] << "\n";
	std::cout << A[WHITE] << "\n";
	std::cout << A[BLACK] << "\n";
#ifdef DEBUG
	for (int c = 0; c < (1 << 3); c++) {
		std::cout << "\"p" << c << "\": [\n";
		//std::cout << T[c].size() << "\n";
		for (const Polygon& t : T[c]) {
			std::cout << "  [\n";
			for (const Pos& p : t) std::cout << "  " << p << ",\n";
			std::cout << "  ],\n";
		}
		std::cout << "],\n";
	}
#endif
	return;
}
//#define TIME
void solve() {
#ifdef TIME
	auto start = std::chrono::high_resolution_clock::now();
#endif
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	//std::cout << std::scientific;
	std::cout.precision(13);
	//freopen("../tests/candle/22.in", "r", stdin);
	//freopen("../tests/candle/22_ans2.txt", "w", stdout);
	std::cin >> N;
	P[0].resize(N);
	for (Pos& p : P[0]) std::cin >> p;
	norm(P[0]);
#ifdef DEBUG
	std::cout << "area(P[0]):: " << area(P[0]) << "\n";
#endif
	std::cin >> K;
	for (int i = 1; i <= K; i++) {
		std::cin >> M;
		P[i].resize(M);
		for (Pos& p : P[i]) std::cin >> p;
#ifdef DEBUG
		std::cout << "area(P[" << i << "]:: " << area(P[i]) << "\n";
#endif
		norm(P[i]);
	}
#ifdef DEBUG
	std::cout << "area(N):: " << area(P[0]) << "\n";
	T[0].push_back(P[0]);
#endif
	std::cin >> Q;
	for (int q = 1; q <= Q; q++) query(q);
#ifdef TIME
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Start time: " << std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count() << " us\n";
	std::cout << "End time: " << std::chrono::duration_cast<std::chrono::microseconds>(end.time_since_epoch()).count() << " us\n";
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Execution time: " << duration.count() << " microseconds\n";
#endif
	return;
}
int main() { solve(); return 0; }//kitpc? 14? candle & candle & candle & shadow

/*
Case #1:
R: 0.0002162129385
G: 0.0000000000000
B: 0.0000000000000
Y: 0.0002898416424
M: 0.0000929948874
C: 0.0000773221836
W: 181929.4966357419326
L: 18070.5027862230490

VS:: 23
f:: 1 se::   (-1000.0000000000000, 989.9496728736790)   (-1000.0000000000000, 900.0000000000000)
f:: 1 se::   (-1000.0000000000000, 900.0000000000000)   (1000.0000000000000, -1000.0000000000000)
f:: 1 se::   (-990.0000000000000, 980.0000000000000)   (-1000.0000000000000, 989.9496728736789)

f:: 0 se::   (-1000.0000000000000, 989.9496981891348)   (-1000.0000000000000, 989.9496728736790)
f:: 0 se::   (-1000.0000000000000, 989.9497234791353)   (-1000.0000000000000, 989.9496981891349)
f:: 0 se::   (-990.0000000000000, 980.0000000000000)   (-1000.0000000000000, 989.9496981891348)
f:: 0 se::   (-990.0000000000000, 980.0000000000000)   (-1000.0000000000000, 989.9497234791353)
f:: 0 se::   (-989.9497234791353, 1000.0000000000000)   (-980.0000000568765, 990.0000000571639)
f:: 0 se::   (-989.9496981891348, 1000.0000000000000)   (-989.9497234791353, 1000.0000000000000)
f:: 0 se::   (-989.9496981891348, 1000.0000000000000)   (-980.0000000395776, 990.0000000397777)
f:: 0 se::   (-989.9496728736789, 1000.0000000000000)   (-989.9496981891348, 1000.0000000000000)
f:: 1 se::   (-989.9496728736789, 1000.0000000000000)   (-980.0000000899036, 990.0000000903583)
f:: 1 se::   (-980.0000000899036, 990.0000000903583)   (-980.0000000366806, 990.0000000368661)
f:: 1 se::   (-980.0000000568765, 990.0000000571639)   (-980.0000000025199, 990.0000000025326)
f:: 1 se::   (-980.0000000395776, 990.0000000397777)   (-980.0000000158310, 990.0000000159110)
f:: 1 se::   (-980.0000000366806, 990.0000000368661)   (-980.0000000118526, 990.0000000119126)
f:: 1 se::   (-980.0000000158310, 990.0000000159110)   (-980.0000000000000, 990.0000000000000)
f:: 1 se::   (-980.0000000118526, 990.0000000119126)   (-980.0000000000000, 990.0000000000000)

f:: 1 se::   (-980.0000000000000, 990.0000000000000)   (501.0000000000000, -499.0000000000000)
f:: 1 se::   (-900.0000000000000, 1000.0000000000000)   (-989.9496728736789, 1000.0000000000000)
f:: 1 se::   (499.0000000000000, -501.0000000000000)   (-990.0000000000000, 980.0000000000000)
f:: 1 se::   (501.0000000000000, -499.0000000000000)   (499.0000000000000, -501.0000000000000)
f:: 1 se::   (1000.0000000000000, -1000.0000000000000)   (-900.0000000000000, 1000.0000000000000)


VS:: 13
f:: 1 se::   (-1000.0000000000000, 989.9496728736790)   (-1000.0000000000000, 900.0000000000000)
f:: 1 se::   (-1000.0000000000000, 900.0000000000000)   (1000.0000000000000, -1000.0000000000000)
f:: 1 se::   (-990.0000000000000, 980.0000000000000)   (-1000.0000000000000, 989.9496728736789)

f:: 0 se::   (-1000.0000000000000, 989.9496981891348)   (-1000.0000000000000, 989.9496728736789)
f:: 0 se::   (-990.0000000000000, 980.0000000000000)   (-1000.0000000000000, 989.9496981891348)
f:: 0 se::   (-989.9496981891348, 1000.0000000000000)   (-980.0000000000000, 990.0000000000000)
f:: 0 se::   (-989.9496728736789, 1000.0000000000000)   (-989.9496981891348, 1000.0000000000000)

f:: 1 se::   (-989.9496728736789, 1000.0000000000000)   (-980.0000000000000, 990.0000000000000)

f:: 1 se::   (-980.0000000000000, 990.0000000000000)   (501.0000000000000, -499.0000000000000)
f:: 1 se::   (-900.0000000000000, 1000.0000000000000)   (-989.9496728736789, 1000.0000000000000)
f:: 1 se::   (499.0000000000000, -501.0000000000000)   (-990.0000000000000, 980.0000000000000)
f:: 1 se::   (501.0000000000000, -499.0000000000000)   (499.0000000000000, -501.0000000000000)
f:: 1 se::   (1000.0000000000000, -1000.0000000000000)   (-900.0000000000000, 1000.0000000000000)


Case #1:
R: 0.0002529000631
G: 0.0000000000000
B: 0.0000000000000
Y: 0.0002531545178
M: 0.0000000000000
C: 0.0000000000000
W: 181929.4967287368199
L: 18070.5027652086283


*/