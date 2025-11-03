#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <deque>
#include <map>
#include <set>
#include <random>
#include <array>
#include <tuple>
#include <complex>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::pair<int, int> pi;
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
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define __FUCK__ ;
#define WHAT_THE_FUCK
//#define DEBUG

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2
#define SEG 3
#define POS 4

#define STRONG 0
#define WEAK 1

int R, N, M;
struct Pos {
	ld x, y;
	int i;
	bool rv;
	Pos(ld X = 0, ld Y = 0) : x(X), y(Y) { i = -1, rv = 0; }
	//bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	//bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	//bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return *this < p || *this == p; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& scalar) const { return { x * scalar, y * scalar }; }
	Pos operator / (const ld& scalar) const { return { x / scalar, y / scalar }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& scale) { x *= scale; y *= scale; return *this; }
	Pos& operator /= (const ld& scale) { x /= scale; y /= scale; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pos rot(ld the) { return { x * cos(the) - y * sin(the), x * sin(the) + y * cos(the) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2l(y, x); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "(" << p.x << ", " << p.y << ")"; return os; }
} p0, p1, key, vec; const Pos O = Pos(0, 0); const Pos INVAL = Pos(INF, INF);
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = cross(d1, d2, d3); return sign(ret); }
bool cmpr(const Pos& p, const Pos& q) {
	//bool f1 = O < p;
	//bool f2 = O < q;
	//if (f1 != f2) return f1;
	//int tq = ccw(O, p, q);
	//return !tq ? p.rv > q.rv : tq > 0;
	ld tp = p.rad();
	ld tq = q.rad();
	return tp == tq ? p.rv > q.rv : tp < tq;
}
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = dot(d1, d3, d2); return !ccw(d1, d2, d3) && sign(ret) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { ld ret = dot(d1, d3, d2); return !ccw(d1, d2, d3) && sign(ret) > 0; }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	//return f1 && f2;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
int inner_check(const Polygon& H, const Pos& p) {//concave
	int sz = H.size(), cnt = 0;
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (cur.y == nxt.y) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y <= p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return (cnt & 1) * 2;
}
ld area(const Polygon& H) {
	ld ret = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) ret += H[i] / H[(i + 1) % sz];
	return ret;
}
struct Seg {
	Pos a, b;
	int i;
	Seg(Pos A = Pos(), Pos B = Pos()) : a(A), b(B) { i = -1; }
	Pos inx(const Seg& o) const { return intersection(a, b, o.a, o.b); }
	Pos p(const ld& rt = .5) const { return a + (b - a) * rt; }
} seg[LEN], frag[LEN * LEN * 10];
typedef std::vector<Seg> Vseg;
Vseg H[100];
Polygon INX[LEN];
void inx_sort(Polygon& INX, const Pos& a) {
	std::sort(INX.begin(), INX.end(), [&](const Pos& p, const Pos& q) -> bool {
		return (a - p).Euc() < (a - q).Euc();
		});
	INX.erase(unique(INX.begin(), INX.end()), INX.end());
	return;
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
	ld area(const ld& lo = 0, const ld& hi = 2 * PI) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		//if (hi < lo) { return green(lo, 2 * PI) + green(0, hi); }
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
		auto the = [&](ld rt) { return norm(q.rad(s + (e - s) * rt)); };
		if (-TOL < hi && hi < 1 + TOL) ret.push_back(the(hi));
		if (zero(det)) return ret;
		if (-TOL < lo && lo < 1 + TOL) ret.push_back(the(lo));
	}
	return ret;
}
int I, I0;
std::map<Pos, Polygon> map_pos;
ld A[LEN * LEN + 10];
Polygon cell[LEN * LEN + 10]; int ci;
std::set<int> cell_i[LEN * LEN + 10];
int P[LEN * LEN + 10];//disjoint set
int find(int i) { return P[i] < 0 ? i : P[i] = find(P[i]); }
bool join(int i, int j) {
	i = find(i), j = find(j);
	if (i == j) return 0;
	if (P[i] < P[j]) P[i] += P[j], P[j] = i;
	else P[j] += P[i], P[i] = j;
	return 1;
}
int V[LEN * LEN * 10];
Vint GS[LEN * LEN * 10];
void dfs(const int& i, int v) {
	V[v] = 1;
	cell[i].push_back(frag[v].a);
	cell_i[i].insert(v);
	for (const int& w : GS[v]) {
		if (V[w]) continue;
		dfs(i, w);
	}
	return;
}
ld par(const Vseg& h) {
	int sz = h.size();
	ld r = 0;
	for (int i = 0; i < sz; i++) {

	}
	return r;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> R >> N;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	for (int i = 0; i < N; i++) P[i].i = i;
	Circle c = Circle(Pos(0, 0), R);
	bool f0 = 1;
	Polygon INXS;
	Vld X;
	for (int i = 0, j; i < N; i++) {
		j = (i + 1) % N;
		const Pos& s = P[i], & e = P[j];
		seg[i].a = s;
		seg[i].b = e;
		seg[i].i = LINE;
		Vld inxs = circle_line_intersections(c, s, e, LINE);
		for (const ld& x : inxs) {
			const Pos p = seg[i].p(x);
			INX[i].push_back(p);
			INX[j].push_back(p);
			INXS.push_back(p);
			X.push_back(x);
		}
	}
	if (f0) {
		ld r = N * 2 * PI;
		for (int i = 0; i < N; i++) {
			const Pos& s = P[i], & e = P[(i + 1) % N];
			r += (s - e).mag();
		}
		std::cout << r << "\n";
		return;
	}
	std::sort(X.begin(), X.end());
	int sz = X.size();
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		seg[i + N].a = c.p(i);
		seg[i + N].b = c.p(j);
		seg[i + N].i = CIRCLE;
	}
	std::sort(INXS.begin(), INXS.end());
	INXS.erase(unique(INXS.begin(), INXS.end()), INXS.end());
	I = 0;
	for (int i = 0; i < M; i++) {
		inx_sort(INX[i], seg[i].a);
		Polygon& v = INX[i];
		int sz = v.size();
		for (int j = 0; j < sz - 1; j++) {
			frag[I] = Seg(v[j], v[j + 1]);
			frag[I].i = I;
			I++;
		}
	}
	I0 = I;
	for (int i = 0; i < M; i++) {
		Polygon& v = INX[i];
		int sz = v.size();
		for (int j = 0; j < sz - 1; j++) {
			frag[I] = Seg(v[j + 1], v[j]);
			frag[I].i = I;
			I++;
		}
	}

	for (int i = 0; i < I; i++) {
		key = frag[i].a;
		vec = frag[i].b - frag[i].a;
		vec.i = frag[i].i;
		map_pos[key].push_back(vec);

		key = frag[i].b;
		vec = frag[i].a - frag[i].b;
		vec.i = frag[i].i;
		vec.rv = 1;
		map_pos[key].push_back(vec);
	}
	for (const Pos& key : INXS) {
		Polygon& v = map_pos[key];
		std::sort(v.begin(), v.end(), cmpr);
		int sz = v.size();
#ifdef DEBUG
		std::cout << "DEBUG:: key:: " << key << "\n";
		std::cout << "DEBUG:: sz:: " << sz << "\n";
		for (int k = 0; k < sz; k++) {
			std::cout << "v[" << k << "]:: " << v[k] << " i:: " << v[k].i << " rv:: " << v[k].rv << "\n";
		}
		std::cout << "FUCK::\n";
#endif
		assert(!(sz & 1));
		for (int j = 0; j < sz; j += 2) {
			Pos cur = v[(j - 1 + sz) % sz], nxt = v[j];
			assert(cur.rv != nxt.rv);
			GS[nxt.i].push_back(cur.i);
		}
	}
	memset(V, 0, sizeof V);
	ci = 0;
	ld AA = 0;
	ld r = 0;
	for (int i = 0; i < I; i++) {
		if (!V[i]) {
			dfs(ci, i);
			A[ci] = area(cell[ci]);
#ifdef DEBUG
			std::cout << "FUCK:: i:: " << i << " sz:: " << cell[ci].size() << "\n";
			//for (const Pos& p : cell[ci]) std::cout << p << "\n";
			std::cout << "FUCK:: A[" << ci << "]:: " << A[ci] << "\n";
#endif
			if (0 == A[ci]) {
				cell[ci].clear();
				cell_i[ci].clear();
				A[ci] = 0;
				ci--;
			}
			ci++;
		}
	}
	return;
}
int main() { solve(); return 0; }//boj
