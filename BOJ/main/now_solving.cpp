#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
typedef long long ll;
typedef double ld;
const ld PI = acos(-1);
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }

int T, N, M;
struct Pos {
	int x, y, i;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return sqrtl(Euc()); }
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
	std::cout.precision(13);
	std::cin >> T; while (T--) {
		ld R = 0, r = 0;
		std::cin >> N; Polygon P(N);
		for (Pos& p : P) std::cin >> p;
		std::cin >> M; Polygon S(M);
		for (int i = 0; i < M; i++) {
			Pos s;
			int k, d; std::cin >> k >> d;
			if (!k) { R += d; continue; }
			const Pos& u = P[k], & v = P[(k + 1) % N];
			s = u;
			if (u.x > v.x) s.x -= d;
			else if (u.x < v.x) s.x += d;
			else if (u.y > v.y) s.y -= d;
			else if (u.y < v.y) s.y += d;
			s.i = 0;
			S.push_back(s);
		}
		M = S.size();
		assert(N >= 4);
		if (P[0].x != P[1].x) { for (Pos& p : P) p = ~p; for (Pos& s : S) s = ~s; }
		if (P[0].x > P[1].x) { for (Pos& p : P) p = -p; for (Pos& s : S) s = -s; }
		Polygon Q;
		for (int i = 0; i < N; i++) if (i != 1) Q.push_back(P[i]);
		for (int i = 0; i < M; i++) Q.push_back(S[i]);
		Polygon H;
		for (const Pos& q : Q) {
			while (H.size() > 1 && ccw(H[H.size() - 2], H.back(), q) <= 0) {
				r -= (H[H.size() - 2] - H.back()).mag();
				H.pop_back();
			}
			r += (q - H.back()).mag();
			H.push_back(q);
			if (!q.i) {
				R += r;
				H.pop_back();
				r -= (q - H.back()).mag();
			}
		}
		std::cout << R << "\n";
	}
	return;
}
int main() { solve(); return 0; }