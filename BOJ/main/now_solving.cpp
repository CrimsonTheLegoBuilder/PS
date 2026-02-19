#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <random>
#include <array>
#include <set>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 1e5 + 1;
const ld TOL = 1e-7;
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }
inline ll sq(const ll& x) { return x * x; }
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

int N, M, T, Q;
struct Pos {
	ll x, y, i;
	Pos(ll x_ = 0, ll y_ = 0, ll i_ = 0) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y, i }; }
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
	Pos unit() const {
		ll g = gcd(std::abs(x), std::abs(y));
		return { x / g, y / g, i };
	}
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
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
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	ll ans = 0;
	for (int i = 0; i < N; i++) {
		Polygon C;
		for (int j = 0; j < N; j++) {
			if (i == j) continue;
			C.push_back((P[j] - P[i]).unit());
		}
		std::sort(C.begin(), C.end(), cmpt);
		int sz = C.size();
		Polygon Q;
		std::set<Pos> S;
		for (int j = 0, k; j < sz; j++) {
			k = j;
			while (k < sz && C[j] == C[k]) k++;
			Pos v = C[j];
			v.i = k - j;
			if (v.i == 0) v.i = 1;
			Q.push_back(v);
			S.insert(v);
			j = k - 1;
		}
		for (const Pos& q : Q) {
			Pos v = ~q;
			auto it = S.find(v);
			if (it != S.end()) {
				ans += q.i * it->i;
			}
		}
	}
	std::cout << ans << "\n";
	return;
}
int main() { solve(); return 0; }//boj3008

