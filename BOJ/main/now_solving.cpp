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
const int LEN = 1e5;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }

#define STRONG 0
#define WEAK 1

int N, M, A;
struct Pos {
	int x, y, i;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} s, e;
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ll area(Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
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
ll f(const Polygon& P, const int& m) {
	Polygon C;
	int sz = P.size();
	for (int i = 0; i < sz; i++) if (P[i].i <= m) C.push_back(P[i]);
	Polygon H = monotone_chain(C);
	return area(H);
}
int bi_search(const Polygon& P, const int& n, const ll& A) {
	int s = 0, e = n - 1;
	if (f(P, e) < A) return -1;
	while (s < e) {
		int m = s + e >> 1;
		ll a = f(P, m);
		if (a == A) return m;
		if (a > A) e = m - 1;
		else s = m + 1;
	}
	return s;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	ll A; std::cin >> N >> A; A <<= 1;
	Polygon P(3), Q(N);
	for (Pos& p : P) std::cin >> p;
	for (Pos& q : Q) std::cin >> q;
	for (int i = 0; i < N; i++) Q[i].i = i;
	for (const Pos& q : Q) P.push_back(q);
	std::sort(P.begin(), P.end());
	int c = bi_search(P, N, A);
	if (c == -1) { std::cout << "draw\n"; return; }
	std::cout << (c & 1 ? "wider\n" : "wapas\n");
	return;
}
int main() { solve(); return 0; }