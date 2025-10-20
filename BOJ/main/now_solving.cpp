#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
typedef long long ll;
typedef double ld;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }

int N;
struct Pos {
	int x, y, i;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
};
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
Polygon graham_scan(Polygon C) {
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
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	while (1) {
		std::cin >> N;
		if (!N) break;
		Polygon P(N);
		for (int i = 0; i < N; i++) std::cin >> P[i], P[i].i = i;
		Polygon H = graham_scan(P);
		int sz = H.size();
		bool f = 1;
		for (int i = 0, s, e; i < sz; i++) {
			s = H[i].i, e = H[(i + 1) % sz].i;
			if ((s + 1) % N == e) continue;
			if (ccw(P[s], P[e], P[(s + 1) % N])) { f = 0; break; }
			if (ccw(P[s], P[e], P[(e - 1 + N) % N])) { f = 0; break; }
		}
		std::cout << (f ? "Yes\n" : "No\n");
	}
	return;
}
int main() { solve(); return 0; }