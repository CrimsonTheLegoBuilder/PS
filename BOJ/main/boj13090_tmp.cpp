#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
typedef long long ll;
typedef long double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-13;
const ld PI = acos(-1);
const int LEN = 30;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }

#define LO x
#define HI y

#define STRONG 0
#define WEAK 1

#define CNT 25

int N, M, K;
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
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y;return is; }
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
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
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
typedef std::vector<Seg> Segs;
ld intersection(const Seg& s1, const Seg& s2, const int& f = 0) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == 2) return a1;
	if (f == 1) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
bool connectable(const Pos& u, const Pos& v, const Polygon& H, const bool& f = 1) {
	int sz = H.size();
	Pos m = (u + v) * .5;
	for (int i = 0; i < sz; i++) {
		const Pos& p0 = H[i], p1 = H[(i + 1) % sz];
		if (intersect(p0, p1, u, v, WEAK)) return 0;
		if (on_seg_weak(u, v, p0)) return 0;
	}
	if (f == 1 && inner_check(H, m) == 2) return 0;
	if (f == 0 && inner_check(H, m) == 0) return 0;
	return 1;
}
struct Event {
	Seg w, e;
	Event(Seg w_ = Seg(), Seg e_ = Seg()) : w(w_), e(e_) {
		if (w.s != w.e && e.s != e.e) {
			Pos v = ~w.dir;
			ld wlo = intersection(w, Seg(e.e, e.e + v), 2);
			ld whi = intersection(w, Seg(e.s, e.s + v), 2);
			ld elo = intersection(e, Seg(w.e, w.e + v), 2);
			ld ehi = intersection(e, Seg(w.s, w.s + v), 2);
			w = Seg(w.p(fit(wlo, 0, 1)), w.p(fit(whi, 0, 1)));
			e = Seg(e.p(fit(ehi, 0, 1)), e.p(fit(elo, 0, 1)));
		}
	}
	bool val() const { return w.s != w.e && e.s != e.e; }
};
typedef std::vector<Event> Events;
struct Info {
	int i;
	ld c;
	Info(int i_ = 0, ld c_ = 0) : i(i_), c(c_) {}
	bool operator < (const Info& x) const { return zero(c - x.c) ? i < x.i : c > x.c; }
};
ld C[LEN * LEN]; int vp;
std::vector<Info> G[LEN * LEN];
std::vector<Info> GW[LEN], GE[LEN];
ld dijkstra(const std::vector<Info> G[], const int& v, const int& g, const int& sz) {
	for (int i = 0; i < sz; i++) C[i] = INF;
	std::priority_queue<Info> Q; Q.push(Info(v, 0));
	C[v] = 0;
	while (Q.size()) {
		Info p = Q.top(); Q.pop();
		if (p.i == g) return C[g];
		if (p.c > C[p.i]) continue;
		for (const Info& w : G[p.i]) {
			ld cost = p.c + w.c;
			if (C[w.i] > cost) {
				C[w.i] = cost;
				Q.push(Info(w.i, cost));
			}
		}
	}
	return C[g];
}
ld dijkstra(
	const Pos& s, const Pos& t,
	const Event& we, const ld& m,
	const Polygon& W, const Polygon& E,
	const Polygon& W_, const Polygon& E_) {
	GW[1].clear(); GE[1].clear();
	Pos w = we.w.p(m), e = we.e.p(m);
	int sz;
	sz = W.size();
	if (connectable(w, s, W_, 0)) GW[1].push_back(Info(0, (s - w).mag()));
	for (int i = 0; i < sz; i++)
		if (connectable(w, W[i], W_, 0))
			GW[1].push_back(Info(i + 2, (w - W[i]).mag()));
	sz = E.size();
	if (connectable(e, t, E_, 0)) GE[1].push_back(Info(0, (t - e).mag()));
	for (int i = 0; i < sz; i++)
		if (connectable(e, E[i], E_, 0))
			GE[1].push_back(Info(i + 2, (e - E[i]).mag()));
	ld dw = dijkstra(GW, 1, 0, N + 5);
	ld de = dijkstra(GE, 1, 0, M + 5);
	return dw + de;
}
ld ternary_search(
	const Pos& s, const Pos& t,
	const Event& we,
	const Polygon& W, const Polygon& E,
	const Polygon& W_, const Polygon& E_) {
	ld lo = 0, hi = 1;
	ld m1, m2, d1, d2;
	int c = CNT;
	while (c--) {
		m1 = (lo + lo + hi) / 3;
		m2 = (lo + hi + hi) / 3;
		d1 = dijkstra(s, t, we, m1, W, E, W_, E_);
		d2 = dijkstra(s, t, we, m2, W, E, W_, E_);
		if (d1 > d2) lo = m1;
		else hi = m2;
	}
	return (d1 + d2) * .5;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	Pos s, t; std::cin >> s >> t;
	std::cin >> N; Polygon W(N); for (Pos& p : W) std::cin >> p;
	std::cin >> M; Polygon E(M); for (Pos& p : E) std::cin >> p;
	if (!eq(W.back().y, 1000)) std::reverse(W.begin(), W.end());
	if (!eq(E.back().y, 0)) std::reverse(E.begin(), E.end());
	Polygon W_ = W; W_.push_back(Pos(0, 1000)); W_.push_back(Pos(0, 0));
	Polygon E_ = E; E_.push_back(Pos(1000, 0)); E_.push_back(Pos(1000, 1000));
	ld b = INF, l = INF, d = -1; Segs B; Events V;
	for (int i = 0; i < N; i++) {
		if (connectable(s, W[i], W_, 0)) {
			d = (s - W[i]).mag();
			GW[0].push_back(Info(i + 2, d));
			GW[i + 2].push_back(Info(0, d));
			G[0].push_back(Info(i + 2, d));
			G[i + 2].push_back(Info(0, d));
		}
		for (int j = i + 1; j < N; j++) {
			if (connectable(W[i], W[j], W_, 0)) {
				d = (W[i] - W[j]).mag();
				GW[i + 2].push_back(Info(j + 2, d));
				GW[j + 2].push_back(Info(i + 2, d));
				G[i + 2].push_back(Info(j + 2, d));
				G[j + 2].push_back(Info(i + 2, d));
			}
		}
	}
	for (int i = 0; i < M; i++) {
		if (connectable(t, E[i], E_, 0)) {
			d = (t - E[i]).mag();
			GE[0].push_back(Info(i + 2, d));
			GE[i + 2].push_back(Info(0, d));
			G[1].push_back(Info(N + i + 2, d));
			G[N + i + 2].push_back(Info(1, d));
		}
		for (int j = i + 1; j < M; j++) {
			if (connectable(E[i], E[j], E_, 0)) {
				d = (E[i] - E[j]).mag();
				GE[i + 2].push_back(Info(j + 2, d));
				GE[j + 2].push_back(Info(i + 2, d));
				G[N + i + 2].push_back(Info(N + j + 2, d));
				G[N + j + 2].push_back(Info(N + i + 2, d));
			}
		}
	}
	for (int n = 0; n < N; n++) {
		const Pos& w0 = W[n], w1 = W[(n + 1) % N];
		Seg w01 = Seg(w0, w1);
		for (int m = 0; m < M; m++) {
			const Pos& e0 = E[m], e1 = E[(m + 1) % M];
			Seg e01 = Seg(e0, e1);
			d = (w0 - e0).mag();
			if (b >= d) {
				if (eq(b, d)) B.push_back(Seg(w0, e0));
				else if (b > d) {
					b = d;
					B.clear(); V.clear(); B.push_back(Seg(w0, e0));
				}
			}
			if (n < N - 1 && between(w0, w1, e0)) {
				d = std::abs(cross(w0, w1, e0)) / (w0 - w1).mag();
				if (b >= d) {
					Pos v = ~w01.dir;
					ld x = intersection(w01, Seg(e0, e0 + v), 2);
					Pos p = w01.p(x);
					if (eq(b, d)) B.push_back(Seg(p, e0));
					else if (b > d) {
						b = d;
						B.clear(); V.clear(); B.push_back(Seg(p, e0));
					}
				}
			}
			if (m < M - 1 && between(e0, e1, w0)) {
				d = std::abs(cross(e0, e1, w0)) / (e0 - e1).mag();
				if (b >= d) {
					Pos v = ~e01.dir;
					ld x = intersection(e01, Seg(w0, w0 + v), 2);
					Pos p = e01.p(x);
					if (eq(b, d)) B.push_back(Seg(w0, p));
					else if (b > d) {
						b = d;
						B.clear(); V.clear(); B.push_back(Seg(w0, p));
					}
				}
			}
			if (n < N - 1 && m < M - 1 &&
				!ccw(w0, w1, e0, e1) && dot(w0, w1, e0, e1) < 0) {
				d = std::abs(cross(w0, w1, e0)) / (w0 - w1).mag();
				if (b >= d) {
					Event ev = Event(w01, e01);
					if (!ev.val()) continue;
					if (eq(b, d)) V.push_back(ev);
					else if (b > d) {
						b = d;
						B.clear(); V.clear(); V.push_back(ev);
					}
				}
			}
		}
	}
	if (B.size()) {
		int sz = B.size();
		int iw = N + M + 2, ie = iw + 1;
		for (int k = 0; k < sz; k++) {
			Pos w = B[k].s, e = B[k].e;
			G[iw].push_back(Info(ie, 0));
			G[ie].push_back(Info(iw, 0));
			if (connectable(s, w, W_, 0)) {
				d = (s - w).mag();
				G[0].push_back(Info(iw, d));
			}
			for (int i = 0; i < N; i++) {
				if (connectable(W[i], w, W_, 0)) {
					d = (W[i] - w).mag();
					G[i + 2].push_back(Info(iw, d));
					G[iw].push_back(Info(i + 2, d));
				}
			}
			if (connectable(t, e, E_, 0)) {
				d = (t - e).mag();
				G[1].push_back(Info(ie, d));
				G[ie].push_back(Info(1, d));
			}
			for (int j = 0; j < M; j++) {
				if (connectable(E[j], e, E_, 0)) {
					d = (E[j] - e).mag();
					G[N + j + 2].push_back(Info(ie, d));
					G[ie].push_back(Info(N + j + 2, d));
				}
			}
			iw += 2; ie += 2;
		}
		d = dijkstra(G, 0, 1, LEN * LEN);
		l = std::min(l, d);
	}
	if (V.size()) {
		for (const Event& we : V) {
			d = ternary_search(s, t, we, W, E, W_, E_);
			l = std::min(l, d);
		}
	}
	std::cout << b << " " << b + l << "\n";
	return;
}
int main() { solve(); return 0; }//boj13090