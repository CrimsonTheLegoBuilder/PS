#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
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
const int LEN = 1e5 + 1;
const ld TOL = 1e-7;
const ll MOD = 1e9 + 7;
const ld PI = acos(-1);
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
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
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
//ll pow_fuck(ll a, ll b) {
//	ll ret = 1;
//	while (b) {
//		if (b & 1) ret = ret * a % MOD;
//		a = a * a % MOD;
//		b >>= 1;
//	}
//	return ret;
//}
//ll powmod(ll a, ll b) {
//	ll res = 1; a %= MOD;
//	assert(b >= 0);
//	for (; b; b >>= 1) {
//		if (b & 1) res = res * a % MOD;
//		a = a * a % MOD;
//	}
//	return res;
//}
struct Info { ll area, l, r; };

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

int N, M, Q;
ll A[LEN], T;
struct Pos {
	int x, y;
	//ll x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	//Pos(ll x_ = 0, ll y_ = 0) : x(x_), y(y_) {}
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
	Pos& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	int Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool cmpt(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ll tq = p / q;
	return !tq ? p.Euc() < q.Euc() : tq > 0;
}
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
void get_area_memo(Pos H[], ll memo[], const int& sz) {
	memo[0] = 0;
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		memo[i + 1] = cross(Pos(0, 0), cur, nxt) + memo[i];//memo[sz] == convex hull's area
	}
	return;
}
void get_round_memo(Polygon& H, ld memo[]) {
	int sz = H.size();
	memo[0] = .0;
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		memo[i + 1] = (cur - nxt).mag() + memo[i];//memo[sz] == convex hull's round
	}
	return;
}
ll area(Pos H[], const int& sz) {
	ll a = 0;
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
ll area(Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) { if (area(H) < 0) std::reverse(H.begin(), H.end()); }
bool inner_check(Pos p0, Pos p1, Pos p2, const Pos& t) {
	if (ccw(p0, p1, p2) < 0) std::swap(p1, p2);
	return ccw(p0, p1, t) >= 0 && ccw(p1, p2, t) >= 0 && ccw(p2, p0, t) >= 0;
}
struct Pdd {
	ld x, y;
	Pdd(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	Pdd operator + (const Pdd& p) const { return { x + p.x, y + p.y }; }
	Pdd operator - (const Pdd& p) const { return { x - p.x, y - p.y }; }
	Pdd operator * (const ld& n) const { return { x * n, y * n }; }
	Pdd operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pdd& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pdd& p) const { return x * p.y - y * p.x; }
	Pdd operator ^ (const Pdd& p) const { return { x * p.x, y * p.y }; }
	Pdd& operator += (const Pdd& p) { x += p.x; y += p.y; return *this; }
	Pdd& operator -= (const Pdd& p) { x -= p.x; y -= p.y; return *this; }
	Pdd& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pdd& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pdd operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pdd rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pdd unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pdd& p1, const Pdd& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pdd& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pdd& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = { 0, 0 };
typedef std::vector<Pdd> Vpdd;
Pdd conv(const Pos& p) { return Pdd(p.x, p.y); }
ld cross(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return sign(cross(d1, d2, d3, d4)); }
ld dot(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pdd& d1, const Pdd& d2, const Pdd& d3, const Pdd& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
Pdd intersection(const Pdd& p1, const Pdd& p2, const Pdd& q1, const Pdd& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld rad(const Pdd& d1, const Pdd& d2, const Pdd& d3) { return rad(d2 - d1, d3 - d1); }
Pdd centroid(const Vpdd& H) {
	Pdd cen = Pdd(0, 0);
	ld a = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		ld a0 = H[i] / H[(i + 1) % sz];
		cen += (H[i] + H[(i + 1) % sz]) * a0;
		a += a0;
	}
	a *= .5;
	cen /= 6;
	if (!zero(a)) cen /= a;
	return cen;
}
ld ternary_search(const Pdd& c, const Pdd& I0, const Pdd& I1, const Pdd& J0, const Pdd& J1, const ld& t = 0) {
	ld l = INF;
	int cnt = 50;
	ld s = 0, e = rad(c, J0, J1);
	Pdd v = J0 - c;
	auto dist = [&](const ld& t) -> ld {
		Pdd v1 = v.rot(norm(s + t));
		Pdd I_ = intersection(c, v1, I0, I1);
		Pdd J_ = intersection(c, v1, J0, J1);
		return (I_ - J_).mag();
		};
	while (cnt--) {
		ld t1 = (s + s + e) / 3;
		ld t2 = (s + e + e) / 3;
		ld l1 = dist(t1);
		ld l2 = dist(t2);
		l = std::min({ l, l1, l2 });
		if (l1 < l2) e = t2;
		else s = t1;
	}
	return l;
}
void query() {
	std::cin >> N; Polygon H(N); for (Pos& p : H) std::cin >> p;
	Vpdd P; for (Pos& p : H) P.push_back(conv(p));
	A[0] = 0;
	for (int i = 0; i < N; i++) A[i + 1] = A[i] + H[i] / H[(i + 1) % N];
	Pdd c = centroid(P);
	ld sh = INF, lg = -1;
	auto area_ = [&](const int& i, const int& j) -> ll {
		if (i <= j) return A[j] - A[i] + H[j] / H[i];
		ll a = A[i] - A[j] + H[i] / H[j];
		return A[N] - a;
		};
	for (int i = 0, j = 1; i < N; i++) {
		while (i != (j + 1) % N && (area_(i, (j + 1) % N) << 1) <= A[N])
			j = (j + 1) % N;
		int j1 = (j + 1) % N;

		Pdd I0 = conv(H[i]), J0 = conv(H[j]), J1 = conv(H[j1]);
		Pdd inx = intersection(c, I0, J0, J1);
		ld l0 = (inx - I0).mag();
		sh = std::min(sh, l0);
		lg = std::max(lg, l0);

		ll a0 = area_(i, j) << 1;

		Pdd I1 = conv(H[(i + 1) % N]);
		Pdd I_ = intersection(I0, I1, J1, c);
		Pdd J_ = intersection(J0, J1, I0, c);
		ld l1 = ternary_search(c, I0, I_, J_, J1);
		sh = std::min(sh, l1);
		lg = std::max(lg, l1);

		I1 = conv(H[(i - 1 + N) % N]);
		if (A[N] != a0) {
			I_ = intersection(I1, I0, J0, c);
			ld l2 = ternary_search(c, I_, I0, J0, J_);
			sh = std::min(sh, l2);
			lg = std::max(lg, l2);
		}
		else {
			J1 = J0;
			J0 = conv(H[(j - 1 + N) % N]);
			I_ = intersection(I1, I0, J0, c);
			J_ = intersection(J0, J1, I0, c);
			ld l2 = ternary_search(c, I_, I0, J0, J_);
			sh = std::min(sh, l2);
			lg = std::max(lg, l2);
		}
	}
	std::cout << sh << " " << lg << "\n";
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> Q; while (Q--) query();
	return;
}
int main() { solve(); return 0; }//11465
//boj 27712 10239 22635 29691 31392 16068

/*

1
8
2 2
-2 -2
2 -2
-2 2
3 0
-3 0
0 3
0 -3

1
10
-2 -4
5 2
7 7
10 -2
9 -5
-5 10
1 -5
4 -9
5 1
7 -7

1
4
4 5
-2 4
1 4
-5 -2

1
5
1 1
2 4
1 4
0 4
-1 1

3
4
0 3
1 3
3 1
3 0
4
-4 0
5 3
0 -4
-1 0
5
4 4
5 0
3 3
3 2
-4 2
5
1 1
2 4
1 4
0 4
-1 1

2
4
-1 -1
1 1
1 -1
-1 1
5
-1 -1
1 1
1 -1
-1 1
0 2


*/

#define BOJ
#ifdef BOJ
#else
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	int Q = 130;
	for (int q = 0; q <= Q; q++) {
		std::string Din = "../../tests/aed/input/input";
		std::string Dout = "../../tests/aed/output/output";
		std::string Qin = Din + std::to_string(q) + ".txt";
		std::string Qout = Dout + std::to_string(q) + ".txt";
		freopen(Qin.c_str(), "r", stdin);
		//freopen("../../tests/aed/input/ret.txt", "w", stdout);
		Pos u;
		std::cin >> N; Polygon H(N); std::cin >> u;
		for (Pos& p : H) std::cin >> p; norm(H);
		Vint I, J;
		for (int i = 0, i1; i < N; i++) {
			i1 = (i + 1) % N;
			Pos p0 = H[i], p1 = H[i1];
			bool vis;
			if (!ccw(u, p0, p1)) {
				if (dot(p0, p1, u) < 0) std::swap(p0, p1);
				vis = 1;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					if (intersect(q0, q1, p1, u)) { vis = 0; break; }
				}
			}
			else {
				if (ccw(u, p0, p1) < 0) std::swap(p0, p1);
				Seg s1 = Seg(p0, p1);
				Vrange V;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					Seg s2 = Seg(q0, q1);
					Range r = range(s1, s2, u);
					if (r.s.den != -1) V.push_back(r);
				}
				vis = 0;
				std::sort(V.begin(), V.end());
				V.push_back(Range(Frac(1), Frac(1)));
				int sz = V.size();
				Frac hi = Frac(0);
				for (const Range& r : V) {
					if (hi < r.s) { vis = 1; break; }
					else hi = std::max(hi, r.e);
				}
			}
			if (vis) I.push_back(i + 1);
		}
		//std::cout << I.size() << "\n";
		//for (int i : I) std::cout << i << " ";
		freopen(Qout.c_str(), "r", stdin);
		std::cin >> M;
		J.resize(M); for (int& j : J) std::cin >> j;
		if (I.size() != M) { std::cout << "fuck::\n"; continue; }
		bool f = 1;
		for (int j = 0; j < M; j++) {
			if (I[j] != J[j]) { f = 0; std::cout << "fuck::\n"; break; }
		}
		if (f) std::cout << "good::\n";
		else std::cout << "fuck::\n";
	}
	return;
}
#endif
//int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068