#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstring>
#include <cassert>
#include <vector>
#include <set>
#include <map>
typedef long long ll;
typedef long long int128;
//typedef __int128 int128;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 1e3 + 5;
const ld TOL = 1e-7;
const ll MOD = 1'000'000'007;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline ll sq(int x) { return (ll)x * x; }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

int N, M;
struct Pos {
	int x, y;
	Pos(int X = 0, int Y = 0) : x(X), y(Y) {}
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
	Pos& operator *= (const int& scale) { x *= scale; y *= scale; return *this; }
	Pos& operator /= (const int& scale) { x /= scale; y /= scale; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	int Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	int quad() const { return y > 0 || y == 0 && x >= 0; }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} u, v; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Pos seg[LEN][2];
Pos B[LEN][12];
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
int collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
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
ll area(const Polygon& H) {
	ll a = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) {
	ll A = area(H); assert(A);
	if (A < 0) std::reverse(H.begin(), H.end());
	auto s = std::min_element(H.begin(), H.end());
	std::rotate(H.begin(), s, H.end());
	return;
}
int128 abs_(int128 x) { if (x < 0) x *= -1; return x; }
int128 gcd(int128 a, int128 b) { while (b) { int128 tmp = a % b; a = b; b = tmp; } return a; }
int sign(const int128& x) { return x < 0 ? -1 : x > 0 ? 1 : 0; }
struct Frac {
	int128 x, den;
	Frac(int128 x_ = 0, int128 den_ = 1) : x(x_), den(den_) {}
	//int f;
	bool operator < (const Frac& o) const {
		int s1 = sign(x) * sign(den);
		int s2 = sign(o.x) * sign(o.den);
		if (s1 != s2) return s1 < s2;
		if (!x) {
			assert(!o.x);
			//if (f != o.f) return f < o.f;
			return 0;
		}
		int128 div1 = abs_(x) / abs_(den);
		int128 div2 = abs_(o.x) / abs_(o.den);
		int128 mod1 = abs_(x) % abs_(den);
		int128 mod2 = abs_(o.x) % abs_(o.den);
		if (div1 == div2) {
			int128 n1 = mod1 * o.den;
			int128 n2 = mod2 * den;
			if (n1 == n2) {
				//if (f != o.f) return f < o.f;
				return 0;
			}
			return s1 > 0 ? n1 < n2 : n1 > n2;
		}
		return s1 > 0 ? div1 < div2 : div1 > div2;
	}
	bool operator == (const Frac& o) const { return x == o.x && den == o.den; }
	friend std::ostream& operator << (std::ostream& os, const Frac& p) { os << p.x << " " << p.den; return os; }
} z(0, 1), o(1, 1);
std::vector<Frac> tmp;
Frac intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) {
	ll a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
	int128 x = (int128)((q2 - q1) / (q1 - p1));
	int128 det = (int128)((q2 - q1) / (p2 - p1));
	int128 d = gcd(abs_(x), abs_(det));
	x /= d; det /= d;
	assert(det && d);
	if (x == 0) det = 1;
	else if (sign(x) * sign(det) < 0 && 0 < x) x *= -1, det *= -1;
	else if (sign(x) * sign(det) > 0) x = abs_(x), det = abs_(det);
	return { x, det };
}
bool block(const int& i, const int& j) {
	bool f0 = intersect(seg[i][0], O, seg[j][0], seg[j][1], WEAK);
	bool f1 = intersect(seg[i][1], O, seg[j][0], seg[j][1], WEAK);
	if (f0 || f1) return 1;
	f0 = intersect(seg[i][0], O, seg[j][0], seg[j][1]);
	f1 = intersect(seg[i][1], O, seg[j][0], seg[j][1]);
	if (f0 && f1) return 1;
	return 0;
}
void shadow(const int& i, const int& j, Frac& s, Frac& e) {
	s = z; e = o;
	if (inside(seg[i][1], O, seg[i][0], seg[j][0]))
		s = intersection(seg[i][0], seg[i][1], O, seg[j][0]);
	if (inside(seg[i][1], O, seg[i][0], seg[j][1]))
		e = intersection(seg[i][0], seg[i][1], O, seg[j][1]);
	return;
}
struct Seg {
	Frac s, e;
	bool operator < (const Seg& o) const { return s == o.s ? e < o.e : s < o.s; }
};
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	std::cin >> N;
	for (int i = 0; i < N; i++) {
		int x, y, l; std::cin >> x >> y >> l;
		Polygon box = { Pos(x, y), Pos(x + l, y), Pos(x + l, y + l), Pos(x, y + l) };
		Pos s = box[0], e = box[0];
		for (const Pos& b : box) {
			if (ccw(O, s, b) < 0) s = b;
			if (ccw(O, e, b) > 0) e = b;
		}
		assert(s / e > 0);
		seg[i][0] = s;
		seg[i][1] = e;
	}
	int cnt = 0;
	for (int i = 0; i < N; i++) {
		//std::cout << "sweep[" << i << "]:: \n";
		std::vector<Seg> V = { { z, z } };
		for (int j = 0; j < N; j++) {
			if (j == i) continue;
			if (block(i, j)) {
				//std::cout << "seg[" << j << "] block seg[" << i << "]\n";
				Frac s, e;
				shadow(i, j, s, e);
				V.push_back({ s, e });
			}
		}
		V.push_back({ o, o });
		std::sort(V.begin(), V.end());
		Frac hi = z;
		bool f = 0;
		for (const Seg& v : V) {
			//std::cout << "v.s:: " << v.s << " v.e:: " << v.e << "\n";
			if (hi < v.s) { f = 1; break; }
			hi = std::max(hi, v.e);
		}
		cnt += f;
	}
	std::cout << cnt << "\n";	
	return;
}
int main() { solve(); return 0; }//boj2397