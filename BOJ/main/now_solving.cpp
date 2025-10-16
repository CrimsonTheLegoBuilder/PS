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

int N;
struct Pos {
	int x, y, i;
	Pos(int x_ = 0, int y_ = 0, int i_ = 0) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { 
		assert(std::abs(x) <= 2e9);
		assert(std::abs(y) <= 2e9);
		assert(std::abs(p.x) <= 2e9);
		assert(std::abs(p.y) <= 2e9);
		return (ll)x * p.x + (ll)y * p.y;
	}
	ll operator / (const Pos& p) const {
		assert(std::abs(x) <= 2e9);
		assert(std::abs(y) <= 2e9);
		assert(std::abs(p.x) <= 2e9);
		assert(std::abs(p.y) <= 2e9);
		return (ll)x * p.y - (ll)y * p.x;
	}
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} s, e;
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
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
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	std::cin >> N; Polygon P(N);
	for (int i = 0; i < N; i++) std::cin >> P[i], P[i].i = i + 1;
	std::cin >> s >> e;
	Polygon C = { s, e }, U, L;
	bool colu = 1, coll = 1;
	int cu = 0, cl = 0;
	for (int i = 0; i < N; i++) {
		Pos p = P[i];
		ll tq = cross(s, e, p), fc = dot(s, e, p);
		if (tq > 0 || (!tq && fc > 0)) {
			U.push_back(p);
			cu++;
			if (tq) colu = 0;
		}
		if (tq < 0 || (!tq && fc < 0)) {
			p.i *= -1;
			L.push_back(p);
			cl++;
			if (tq) coll = 0;
		}
		C.push_back(p);
	}
	if (colu && coll) { std::cout << "NO\n"; return; }
	if (!cu || !cl) {
		std::sort(P.begin(), P.end());
		std::cout << "YES\n";
		for (int i = 0; i < N - 1; i++) std::cout << P[i].i << " " << P[i + 1].i << "\n";
		return;
	}
	Polygon H = monotone_chain(C);
	int sz = H.size();
	Pos u, l;
	bool f = 0;
	//for (const Pos& p : H) std::cout << "H:: " << p.i << "\n";
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		if (H[i].i * H[j].i < 0) {
			if (H[i].i > 0) u = H[i], l = H[j];
			else l = H[i], u = H[j];
			f = 1;
			assert(!intersect(u, l, s, e));
			break;
		}
	}
	//std::cout << "fuck::\n";
	if (!f) { std::cout << "NO\n"; return; }
	assert(u.i); assert(l.i);
	assert(u.i * l.i < 0);
	//for (int i = 0; i < cu; i++) assert(U[i] != u);
	//for (int i = 0; i < cl; i++) assert(L[i] != l);
	f = 0;
	for (int i = 0; i < cu; i++) if (U[i] == u) { f = 1; std::swap(U[i], U[0]); break; }
	assert(f);
	std::sort(U.begin() + 1, U.end(), [&](const Pos& p, const Pos& q) -> bool {
		ll tq = ccw(U[0], p, q);
		return !tq ? p.Euc() < q.Euc() : tq > 0;
		});
	f = 0;
	for (int i = 0; i < cl; i++) if (L[i] == l) { f = 1; std::swap(L[i], L[0]); break; }
	assert(f); 
	std::sort(L.begin() + 1, L.end(), [&](const Pos& p, const Pos& q) -> bool {
		ll tq = ccw(L[0], p, q);
		return !tq ? p.Euc() < q.Euc() : tq > 0;
		});
	std::cout << "YES\n";
	std::cout << u.i << " " << -l.i << "\n";
	for (int i = 0; i < cu - 1; i++) std::cout << U[i].i << " " << U[i + 1].i << "\n";
	//for (int i = 0; i < cl - 1; i++) std::cout << -L[i].i << " " << -L[i + 1].i << "\n";
	return;
}
int main() { solve(); return 0; }