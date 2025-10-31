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
inline int sign(ll x) { return x < 0 ? -1 : !!x; }

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
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	std::sort(C.begin(), C.end());
	int sz = C.size();
	if (sz <= 2) return C;
	bool f = 1;
	for (int i = 1; i < sz - 1; i++) if (ccw(C[i - 1], C[i], C[i + 1])) { f = 0; break; }
	if (f) return C;
	for (int i = 0; i < C.size(); i++) {
		while (H.size() > 1 && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) < 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	H.pop_back();
	int s = H.size() + 1;
	for (int i = C.size() - 1; i >= 0; i--) {
		while (H.size() > s && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) < 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	H.pop_back();
	return H;
}
ll f(Polygon P, const Polygon& Q, const int& m) {
	for (int i = 0; i <= m; i++) P.push_back(Q[i]);
	Polygon H = graham_scan(P);
	return area(H);
}
int bi_search(const Polygon& P, const Polygon& Q, const ll& A) {
	int s = 0, e = Q.size();
	if (f(P, Q, e) < A) return -1;
	while (s < e) {
		int m = s + e >> 1;
		ll a = f(P, Q, m);
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
	ll A;
	std::cin >> N >> A; Polygon P(3), Q(N);
	int c = bi_search(P, Q, A);
	if (c == -1) { std::cout << "draw\n"; return; }
	std::cout << (c & 1 ? "wider\n" : "wapas\n");
	return;
}
int main() { solve(); return 0; }