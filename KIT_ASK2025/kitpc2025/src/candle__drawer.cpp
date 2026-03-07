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
const ld TOL = 1e-9;
const int LEN = 105;
const ld PI = acos(-1);
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
#define DEBUG
#ifdef DEBUG
#endif

//#define ASSERT
#ifdef ASSERT
#endif

#define POLYGON_CHECK

#define NAIVE
#ifndef NAIVE
#define FAST
#endif
//=======================DEBUG MACRO=======================//

char CC[8] = {
	'L',
	'R',
	'G',
	'Y',
	'B',
	'M',
	'C',
	'W',
};

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
	Pos rot(const ld& t) { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "  (" << p.x << ", " << p.y << ")"; return os; }
} L[1 << 3]; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon P[LEN];
std::vector<Polygon> T[1 << 3];
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
	Polygon ret = C;
	for (int i = 0; i < sz; i++) {
		Pos b1 = clip[i], b2 = clip[(i + 1) % sz];
		ret = convex_cut(ret, b1, b2);
	}
	return ret;
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
	Pos dir() const { return (s - e).unit(); }
	//bool operator < (const Seg& l) const {
	//	Pos v0 = dir();
	//	Pos v1 = l.dir();
	//	bool f0 = O < v0;
	//	bool f1 = O < v1;
	//	if (f0 != f1) return f0;
	//	if (collinear(s, e, l.s, l.e)) return s < l.s;
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
bool inner_check(const Polygon& H, const Pos& q, const Pos& d = Pos(0, 0)) {
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
Polygon intersection(const Polygon& a, const Polygon& b) {
	Polygon ret = sutherland_hodgman(a, b);
	//ret = graham_scan(ret);
	//if (zero(area(ret))) return {};
	return ret;
}
Polygon intersection(const Polygon& a, const Polygon& b, const Polygon& c) {
	Polygon d = sutherland_hodgman(a, b);
	d = graham_scan(d);
	if (zero(area(d))) return {};
	Polygon ret = sutherland_hodgman(d, c);
	return ret;
	//ret = graham_scan(ret);
	//if (zero(area(ret))) return {};
	//return ret;
}
Vld intersections(const Seg& l, const Polygon& H) {
	int sz = H.size();
	Vld ret;
	//std::cout << "l.s:: " << l.s << " ";
	//std::cout << "l.e:: " << l.e << "\n";
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = H[i], & p1 = H[(i + 1) % sz];
		//std::cout << "p0:: " << p0 << " ";
		//std::cout << "p1:: " << p1 << "\n";
		Seg k = Seg(p0, p1);
		if (collinear(l.s, l.e, p0, p1)) {
			if (dot(l.s, l.e, p0, p1) < 0) return {};
			for (const Pos& p : { p0, p1 }) {
				ld ix = projection(l.s, l.e, p);
				//std::cout << "prj:: " << ix << "\n";
				ret.push_back(fit(ix));
			}
		}
		else if (ccw(l.s, l.e, p0) * ccw(l.s, l.e, p1) <= 0) {
			ld ix = intersection(l, k, WEAK);
			//std::cout << "ix:: " << ix << "\n";
			ret.push_back(fit(ix));
		}
	}
	std::sort(ret.begin(), ret.end());
	ret.erase(unique(ret.begin(), ret.end(), eq), ret.end());
	return ret;
}
//Vld intersections(const Seg& l, const int& c) {
//	Vld ret;
//	for (const Polygon& t : T[c]) {
//		Vld tmp = intersections(l, t);
//		ret.insert(ret.end(), tmp.begin(), tmp.end());
//	}
//	return ret;
//}
//bool inner_check_seg(const int& c, const Seg& l, const ld& x) {
//	Pos q = l.p(x), d = l.dir();
//	Pos pet = ~d;//centripetal
//	for (const Polygon& t : T[c]) {
//		assert(t.size() == 3);
//		if (inner_check(t, q, d)) return 1;
//	}
//	return 0;
//}
struct Event {
	ld x;
	int f;
	bool operator < (const Event& o) const { return eq(x, o.x) ? f < o.f : sign(x - o.x) < 0; }
};
typedef std::vector<Event> Ve;
void query(const int& q) {
	std::cin >> L[RED] >> L[GREEN] >> L[BLUE];
	memset(A, 0, sizeof A);
	memset(I, -1, sizeof I);
	memset(C, 0, sizeof C);
	for (int c = 1; c < (1 << 3); c++) T[c].clear();
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
					Pos s = w.p(hi);
					Pos e = w.p(p.LO);
					Polygon tri = { L[c], s, e };
					T[c].push_back(tri);
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
						Pos s = w.p(hi);
						Pos e = w.p(p.LO);
						Polygon tri = { L[c], s, e };
						T[c].push_back(tri);
						hi = p.HI;
					}
					else hi = std::max(hi, p.HI);
				}
			}
		}
	}
	if (!I[RED] && !I[GREEN] && !I[BLUE]) {//R & G & B
#ifdef NAIVE //27N^3 -> N==200:216,000,000
		for (const Polygon& R : T[RED]) {
			for (const Polygon& G : T[GREEN]) {
				for (const Polygon& B : T[BLUE]) {
					Polygon ix = intersection(R, G, B);
					A[WHITE] += area(ix);
#ifdef POLYGON_CHECK
					if (ix.size()) T[WHITE].push_back(ix);
#endif
				}
			}
		}
#endif
#ifdef FAST //3*3*2N^2log2N^2 -> N==200:11,520,000~
		Segs VS;
		for (int i = 0; i < 3; i++) {
			int c0 = 1 << i;
			int c1 = (1 << ((i + 1) % 3));
			int c2 = (1 << ((i + 2) % 3));
			for (const Polygon& T0 : T[c0]) {
				for (int j = 0; j < 3; j++) {
					const Pos& p0 = T0[j], & p1 = T0[(j + 1) % 3];
					Seg l = Seg(p0, p1);
					Polygon vp;
					Ve ve = { { 0, 1 }, { 1, -1 } };
					Vld vx = { 0, 1 };
					for (const int& c : { c1, c2 }) {
						for (const Polygon& t : T[c]) {
							Vld tmp = intersections(l, t);
#ifdef DEBUG
							std::cout << "tmp.size():: " << tmp.size() << "\ntmp::\n";
							for (const ld& x : tmp) std::cout << x << " ";
							std::cout << "\ntmp\n";
#endif
							if (tmp.size() > 1) {
								assert(tmp.size() == 2);
								ld s = tmp[0];
								ld e = tmp[1];
								ve.push_back({ s, 1 });
								ve.push_back({ e, -1 });
								vx.push_back(s);
								vx.push_back(e);
							}
						}
					}
					std::sort(ve.begin(), ve.end());
					std::sort(vx.begin(), vx.end());
					vx.erase(unique(vx.begin(), vx.end()), vx.end());
					int szr = ve.size(), szx = vx.size(), cnt = 0;
					for (int x = 0, k = 0; x < szx - 1; x++) {
						const ld& s = vx[x], e = vx[x + 1];
						while (k < szr && eq(ve[k].x, s)) { cnt += ve[k].f; k++; }
						if (cnt > 2) vp.push_back(Pos(s, e));
					}
					for (const Pos& se : vp) VS.push_back(Seg(l.p(se.x), l.p(se.y)));
				}
			}
		}
		std::sort(VS.begin(), VS.end());
		VS.erase(unique(VS.begin(), VS.end()), VS.end());
		for (const Seg& se : VS) A[WHITE] += se.green();
#endif
	}
	for (int i = 0; i < 3; i++) {//R & G, G & B, B & R  3*9*N^2 -> N==200:1,080,000
		int c1 = (1 << ((i + 1) % 3));
		int c2 = (1 << ((i + 2) % 3));
		if (!I[c1] && !I[c2]) {
			int c = c1 | c2;
			for (const Polygon& T1 : T[c1]) {
				for (const Polygon& T2 : T[c2]) {
					Polygon ix = intersection(T1, T2);
					A[c] += area(ix);
#ifdef POLYGON_CHECK
					if (ix.size()) T[c].push_back(ix);
#endif
				}
			}
		}
	}
	for (int c = 1; c < (1 << 3); c <<= 1) {//R, G, B  3*3*N -> N==200:1,800
		if (!I[c]) {
			for (const Polygon& T0 : T[c]) A[c] += area(T0);
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
	std::cout << A[RED] << "\n";
	std::cout << A[GREEN] << "\n";
	std::cout << A[BLUE] << "\n";
	std::cout << A[YELLOW] << "\n";
	std::cout << A[MAGENTA] << "\n";
	std::cout << A[CYAN] << "\n";
	std::cout << A[WHITE] << "\n";
	std::cout << A[BLACK] << "\n";
#ifdef DEBUG
	//for (int c = 0; c < (1 << 3); c++) {
	//	std::cout << CC[c] << " = [\n";
	//	//std::cout << "\"p" << c << "\": [\n";
	//	//std::cout << T[c].size() << "\n";
	//	for (const Polygon& t : T[c]) {
	//		std::cout << "  [\n";
	//		for (const Pos& p : t) std::cout << "  " << p << ",\n";
	//		std::cout << "  ],\n";
	//	}
	//	std::cout << "],\n";
	//}
	for (int c = 1; c <= 1; c++) {
		std::cout << CC[c] << " = [\n";
		//std::cout << "\"p" << c << "\": [\n";
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
void solve() {
	auto start = std::chrono::high_resolution_clock::now();
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	//std::cout << std::scientific;
	std::cout.precision(8);
	freopen("../tests/candle/in/80.in", "r", stdin);
	freopen("../tests/candle/test_80.txt", "w", stdout);
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
	//std::cin >> Q;
	Q = 1;
	for (int q = 1; q <= Q; q++) query(q);
	//auto end = std::chrono::high_resolution_clock::now();
	//std::cout << "Start time: " << std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count() << " us\n";
	//std::cout << "End time: " << std::chrono::duration_cast<std::chrono::microseconds>(end.time_since_epoch()).count() << " us\n";
	//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	//std::cout << "Execution time: " << duration.count() << " microseconds\n";
	return;
}
int main() { solve(); return 0; }//kitpc? 14? candle & candle & candle & shadow


/*

4
-4 -4
4 -4
4 4
-4 4
1
4
-1 -1
1 -1
1 1
-1 1
1
3 3 -3 3 -3 -3

100
-966 -1000
-879 -980
-794 -960
-711 -940
-630 -920
-551 -900
-474 -880
-399 -860
-326 -840
-255 -820
-186 -800
-119 -780
-54 -760
9 -740
70 -720
129 -700
186 -680
241 -660
294 -640
345 -620
394 -600
441 -580
486 -560
529 -540
570 -520
609 -500
646 -480
681 -460
714 -440
745 -420
774 -400
801 -380
826 -360
849 -340
870 -320
889 -300
906 -280
921 -260
934 -240
945 -220
955 -200
964 -180
972 -160
979 -140
985 -120
990 -100
994 -80
997 -60
999 -40
1000 -20
1000 20
999 40
997 60
994 80
990 100
985 120
979 140
972 160
964 180
955 200
945 220
934 240
921 260
906 280
889 300
870 320
849 340
826 360
801 380
774 400
745 420
714 440
681 460
646 480
609 500
570 520
529 540
486 560
441 580
394 600
345 620
294 640
241 660
186 680
129 700
70 720
9 740
-54 760
-119 780
-186 800
-255 820
-326 840
-399 860
-474 880
-551 900
-630 920
-711 940
-794 960
-879 980
-966 1000
33
3
100 -97
101 -96
100 -95
3
100 -91
101 -90
100 -89
3
100 -85
101 -84
100 -83
3
100 -79
101 -78
100 -77
3
100 -73
101 -72
100 -71
3
100 -67
101 -66
100 -65
3
100 -61
101 -60
100 -59
3
100 -55
101 -54
100 -53
3
100 -49
101 -48
100 -47
3
100 -43
101 -42
100 -41
3
100 -37
101 -36
100 -35
3
100 -31
101 -30
100 -29
3
100 -25
101 -24
100 -23
3
100 -19
101 -18
100 -17
3
100 -13
101 -12
100 -11
3
100 -7
101 -6
100 -5
3
100 -1
101 0
100 1
3
100 5
101 6
100 7
3
100 11
101 12
100 13
3
100 17
101 18
100 19
3
100 23
101 24
100 25
3
100 29
101 30
100 31
3
100 35
101 36
100 37
3
100 41
101 42
100 43
3
100 47
101 48
100 49
3
100 53
101 54
100 55
3
100 59
101 60
100 61
3
100 65
101 66
100 67
3
100 71
101 72
100 73
3
100 77
101 78
100 79
3
100 83
101 84
100 85
3
100 89
101 90
100 91
3
100 95
101 96
100 97
1
800 -301 800 29 800 359



9888.9209382589906
1775.9727990110405
10016.9697187361307
110738.8909080252051
102497.8939882982522
110610.8421275485307
2396034.2255939878523
1556.2839261339977

9888.92093828
1647.92401855
9888.92093829
110738.89090802
102497.89398828
110738.89090802
2396034.22559398
1684.33270657


100
800 -380 800 -50 800 280
800 -379 800 -49 800 281
800 -378 800 -48 800 282
800 -377 800 -47 800 283
800 -376 800 -46 800 284
800 -375 800 -45 800 285
800 -374 800 -44 800 286
800 -373 800 -43 800 287
800 -372 800 -42 800 288
800 -371 800 -41 800 289
800 -370 800 -40 800 290
800 -369 800 -39 800 291
800 -368 800 -38 800 292
800 -367 800 -37 800 293
800 -366 800 -36 800 294
800 -365 800 -35 800 295
800 -364 800 -34 800 296
800 -363 800 -33 800 297
800 -362 800 -32 800 298
800 -361 800 -31 800 299
800 -360 800 -30 800 300
800 -359 800 -29 800 301
800 -358 800 -28 800 302
800 -357 800 -27 800 303
800 -356 800 -26 800 304
800 -355 800 -25 800 305
800 -354 800 -24 800 306
800 -353 800 -23 800 307
800 -352 800 -22 800 308
800 -351 800 -21 800 309
800 -350 800 -20 800 310
800 -349 800 -19 800 311
800 -348 800 -18 800 312
800 -347 800 -17 800 313
800 -346 800 -16 800 314
800 -345 800 -15 800 315
800 -344 800 -14 800 316
800 -343 800 -13 800 317
800 -342 800 -12 800 318
800 -341 800 -11 800 319
800 -340 800 -10 800 320
800 -339 800 -9 800 321
800 -338 800 -8 800 322
800 -337 800 -7 800 323
800 -336 800 -6 800 324
800 -335 800 -5 800 325
800 -334 800 -4 800 326
800 -333 800 -3 800 327
800 -332 800 -2 800 328
800 -331 800 -1 800 329
800 -330 800 0 800 330
800 -329 800 1 800 331
800 -328 800 2 800 332
800 -327 800 3 800 333
800 -326 800 4 800 334
800 -325 800 5 800 335
800 -324 800 6 800 336
800 -323 800 7 800 337
800 -322 800 8 800 338
800 -321 800 9 800 339
800 -320 800 10 800 340
800 -319 800 11 800 341
800 -318 800 12 800 342
800 -317 800 13 800 343
800 -316 800 14 800 344
800 -315 800 15 800 345
800 -314 800 16 800 346
800 -313 800 17 800 347
800 -312 800 18 800 348
800 -311 800 19 800 349
800 -310 800 20 800 350
800 -309 800 21 800 351
800 -308 800 22 800 352
800 -307 800 23 800 353
800 -306 800 24 800 354
800 -305 800 25 800 355
800 -304 800 26 800 356
800 -303 800 27 800 357
800 -302 800 28 800 358

800 -300 800 30 800 360
800 -299 800 31 800 361
800 -298 800 32 800 362
800 -297 800 33 800 363
800 -296 800 34 800 364
800 -295 800 35 800 365
800 -294 800 36 800 366
800 -293 800 37 800 367
800 -292 800 38 800 368
800 -291 800 39 800 369
800 -290 800 40 800 370
800 -289 800 41 800 371
800 -288 800 42 800 372
800 -287 800 43 800 373
800 -286 800 44 800 374
800 -285 800 45 800 375
800 -284 800 46 800 376
800 -283 800 47 800 377
800 -282 800 48 800 378
800 -281 800 49 800 379

*/