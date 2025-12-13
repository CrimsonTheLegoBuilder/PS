#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
const ll INF = 1e17;
const int LEN = 1e5 + 1;
const ld TOL = 1e-7;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }
inline ll sq(const ll& x) { return x * x; }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

const int N_ = 1 << 16;
const int K_ = 1 << 8 | 1 << 6;

int N, M, K, T, Q;
struct Pos {
	int x, y;
	//ll x, y;
	int i, j;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1, int j_ = -1) : x(x_), y(y_), i(i_), j(j_) {}
	//Pos(ll x_ = 0, ll y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon H[K_];
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpx_rvs(const Pos& p, const Pos& q) { return p.x == q.x ? p.y > q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpy_rvs(const Pos& p, const Pos& q) { return p.y == q.y ? p.x > q.x : p.y < q.y; }
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
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
int collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	//std::sort(C.begin(), C.end());
	//C.erase(unique(C.begin(), C.end()), C.end());
	if (C.size() <= 2) { for (const Pos& pos : C) H.push_back(pos); }
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
int top[K_], bot[K_], lb[K_], rb[K_];//convex hull rotating calipers jaw data
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N >> K;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	std::sort(P.begin(), P.end());
	for (int i = 0; i < N; i++) P[i].i = i;
	Vbool F(N, 0);
	Polygon C = P;
	int k_ = 0;
	for (int k = 0; k < K + 1; k++) {//O(NK)
		Polygon tmp;
		for (const Pos& c : C) if (!F[c.i]) tmp.push_back(c);
		if (tmp.empty()) break;
		k_++;
		H[k] = monotone_chain(tmp);
		for (const Pos& h : H[k]) F[h.i] = 1;
		C.clear();
		for (const Pos& t : tmp) if (!F[t.i]) C.push_back(t);
	}
	std::sort(P.begin(), P.end(), cmpx_rvs);
	for (int i = 0, j; i < K + 1; i++) {
		j = N - i - 1;
		bot[i] = P[i].i;
		top[i] = P[j].i;
	}
	Pos main_vec = Pos(0, -1);
	Pos vbot = Pos(0, -1), vtop = Pos(0, 1);
	while (1) {

	}
}
int main() { solve(); return 0; }//boj26108