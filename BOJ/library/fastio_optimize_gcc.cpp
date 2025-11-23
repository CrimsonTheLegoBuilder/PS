#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#define _CRT_SECURE_NO_WARNINGS
//#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <random>
#include <array>
#include <tuple>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 3e4 + 5;
const ld TOL = 1e-7;
const ll MOD = 1e9 + 7;
const ld PI = acos(-1);
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }

#include <unistd.h>

int iptr = 0, left = 0;
char ibuf[1 << 20];

int read() {
	int x = 0, s = 1;
	bool started = false;
	for (;;) {
		if (iptr >= left) {
			left = read(0, ibuf, sizeof ibuf);
			iptr = 0;
		}
		if (!left) break;
		char c = ibuf[iptr++];
		if (c == '-') {
			s = -1;
		}
		else if (c >= '0' && c <= '9') {
			x = x * 10 + (c - '0');
			started = true;
		}
		else if (started) {
			break;
		}
	}
	return s > 0 ? x : -x;
}

#include <unistd.h>
namespace fast_io {
	const int BUF_SIZE = 1 << 20;
	char ibuf[BUF_SIZE];
	int iptr = 0, left = 0;

	int read_int() {
		int x = 0, s = 1;
		bool started = false;
		for (;;) {
			if (iptr >= left) {
				left = read(0, ibuf, sizeof ibuf);
				iptr = 0;
			}
			if (!left) break;
			char c = ibuf[iptr++];
			if (c == '-') {
				s = -1;
			}
			else if (c >= '0' && c <= '9') {
				x = x * 10 + (c - '0');
				started = true;
			}
			else if (started) {
				break;
			}
		}
		return s > 0 ? x : -x;
	}

	char obuf[BUF_SIZE];
	int optr = 0;

	void flush() {
		if (optr) {
			write(1, obuf, optr);
			optr = 0;
		}
	}

	struct Flusher {
		~Flusher() {
			flush();
		}
	} flusher;

	void pc(char c) {
		if (optr == BUF_SIZE) {
			flush();
		}
		obuf[optr++] = c;
	}

	void print_int(int x) {
		if (x == 0) {
			pc('0');
			return;
		}
		if (x < 0) {
			pc('-');
			x = -x;
		}
		char temp[10];
		int idx = 0;
		while (x) {
			temp[idx++] = x % 10 + '0';
			x /= 10;
		}
		while (idx--) {
			pc(temp[idx]);
		}
	}
}

#define LO x
#define HI y

#define STRONG 0
#define WEAK 1

int N, Q, R;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	[[gnu::const, gnu::hot]]
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	[[gnu::const, gnu::hot]]
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	[[gnu::const, gnu::hot]]
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	[[gnu::const, gnu::hot]]
	bool operator <= (const Pos& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	[[gnu::const, gnu::hot]]
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	[[gnu::const, gnu::hot]]
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	[[gnu::const, gnu::hot]]
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	[[gnu::const, gnu::hot]]
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	[[gnu::const, gnu::hot]]
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	[[gnu::const, gnu::hot]]
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	[[gnu::const, gnu::hot]]
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	[[gnu::const, gnu::hot]]
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	[[gnu::const, gnu::hot]]
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	[[gnu::const, gnu::hot]]
	Pos& operator *= (const int& n) { x *= n; y *= n; return *this; }
	[[gnu::const, gnu::hot]]
	Pos& operator /= (const int& n) { x /= n; y /= n; return *this; }
	[[gnu::const, gnu::hot]]
	Pos operator - () const { return { -x, -y }; }
	[[gnu::const, gnu::hot]]
	Pos operator ~ () const { return { -y, x }; }
	[[gnu::const, gnu::hot]]
	Pos operator ! () const { return { y, x }; }
	[[gnu::const, gnu::hot]]
	ll xy() const { return (ll)x * y; }
	[[gnu::const, gnu::hot]]
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	[[gnu::const, gnu::hot]]
	int Man() const { return std::abs(x) + std::abs(y); }
	[[gnu::const, gnu::hot]]
	ld mag() const { return hypot(x, y); }
	[[gnu::const, gnu::hot]]
	ld rad() const { return atan2(y, x); }
	[[gnu::const, gnu::hot]]
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	//friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	//friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} S[LEN], E[LEN]; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
[[gnu::const, gnu::hot]]
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
[[gnu::const, gnu::hot]]
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
[[gnu::const, gnu::hot]]
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
[[gnu::const, gnu::hot]]
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
[[gnu::const, gnu::hot]]
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
[[gnu::const, gnu::hot]]
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
[[gnu::const, gnu::hot]]
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
[[gnu::const, gnu::hot]]
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
[[gnu::const, gnu::hot]]
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
[[gnu::const, gnu::hot]]
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
bool inside(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q, const int& f = STRONG) {
	if (ccw(p0, p1, p2) < 0) return ccw(p0, p1, q) >= f || ccw(p1, p2, q) >= f;
	return ccw(p0, p1, q) >= f && ccw(p1, p2, q) >= f;
}
[[gnu::const, gnu::hot]]
ld intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) {
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	return ((q2 - q1) / (q1 - p1)) / det;
}
[[gnu::const, gnu::hot]]
bool inner_check(const Polygon& H, const Pos& q, const Pos& d = Pos(0, 0)) {//convex
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos& p1 = H[i], & p2 = H[(i + 1) % sz];
		if (ccw(p1, p2, q) < 0) return 0;
		//if (on_seg_strong(p1, p2, q)) return 2;
	}
	return 1;
}
[[gnu::const, gnu::hot]]
Pos shadow(Pos l, Pos p1, Pos p2, Pos q1, Pos q2, const int& n) {
	if (p1.y == p2.y) { l = ~l, p1 = ~p1, p2 = ~p2, q1 = ~q1, q2 = ~q2; }
	if (p1.y > p2.y) { l = -l, p1 = -p1, p2 = -p2, q1 = -q1, q2 = -q2; }
	int diff = p1.y;
	l.y -= diff, p1.y -= diff, p2.y -= diff, q1.y -= diff, q2.y -= diff;
	//assert(ccw(l, p1, p2) > 0);
	//assert(ccw(l, q1, q2) > 0);
	if (!inside(p2, l, p1, q1, WEAK) && !inside(p2, l, p1, q2, WEAK)) {
		if (intersect(l, p1, q1, q2) && intersect(l, p2, q1, q2)) return Pos(0, n);
		else return INVAL;
	}
	Polygon tri = { p1, p2, l };
	bool in1 = inner_check(tri, q1), in2 = inner_check(tri, q2);
	if (!in1 && !in2) return INVAL;
	ld r1 = 0, r2 = n;
	int lo = 0, hi = n;
	Pos x;
	if (in1 && in2) {
		r1 = intersection(p1, p2, l, q1);
		r1 *= n;
		x = Pos(p1.x, (int)(r1 + TOL));
		if (!ccw(l, q1, x)) lo = x.y;
		else lo = x.y + 1;
		r2 = intersection(p1, p2, l, q2);
		r2 *= n;
		x = Pos(p1.x, (int)(r2 + TOL));
		if (!ccw(l, q2, x)) hi = x.y;
		else hi = x.y;
	}
	else if (in1) {
		r1 = intersection(p1, p2, l, q1);
		r1 *= n;
		x = Pos(p1.x, (int)(r1 + TOL));
		if (!ccw(l, q1, x)) lo = x.y;
		else lo = x.y + 1;
	}
	else if (in2) {
		r2 = intersection(p1, p2, l, q2);
		r2 *= n;
		x = Pos(p1.x, (int)(r2 + TOL));
		if (!ccw(l, q2, x)) hi = x.y;
		else hi = x.y;
	}
	else lo = hi = -1;
	if (lo > hi) lo = hi = -1;
	return Pos(lo, hi);
}
void solve() {
	//std::cin.tie(0)->sync_with_stdio(0);
	//std::cout.tie(0);
	Pos J;
	//std::cin >> N >> R >> J;
	N = read();
	R = read();
	J.x = read();
	J.y = read();
	for (int i = 0; i < R; i++) {
		//int n; std::cin >> n;
		int n; n = read();
		Polygon P(n);
		//for (Pos& p : P) std::cin >> p;
		for (Pos& p : P) p.x = read(), p.y = read();
		Pos s = P[0], e = P[0];
		for (const Pos& p : P) {
			if (ccw(J, s, p) < 0) s = p;
			if (ccw(J, e, p) > 0) e = p;
		}
		S[i] = s; E[i] = e;
	}
	int ret = 0;
	Polygon B = { Pos(0, 0), Pos(N, 0), Pos(N, N), Pos(0, N) };
	for (int b = 0; b < 4; b++) {
		Polygon vp = { Pos(0, 0) };
		Pos s = B[b], e = B[(b + 1) % 4];
		bool f = 1;
		for (int i = 0; i < R; i++) {
			Pos se = shadow(J, s, e, S[i], E[i], N);
			if (se.LO == -1) continue;
			if (se.LO == 0) f = 0;
			//if (inside(E[i], J, S[i], s)) f = 0;
			vp.push_back(se);
		}
		if (f) ret++;
		vp.emplace_back(N, N);
		std::sort(vp.begin(), vp.end());
		int hi = 0;
		for (const Pos& r : vp) {
			if (hi < r.LO) {
				ret += r.LO - hi - 1;
				hi = r.HI;
			}
			else hi = std::max(hi, r.HI);
		}
	}
	//std::cout << ret << "\n";
	printf("%d\n", ret);
	for (i = 0; i < N; i++) {
		fast_io::print_int(P[i].x + h.x);
		fast_io::pc(' ');
		fast_io::print_int(P[i].y + h.y);
		fast_io::pc('\n');
	}
	return;
}
int main() { solve(); return 0; }