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

int N;
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
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	std::cin >> N; Polygon P(N); for (int i = 0; i < N; i++) std::cin >> P[i], P[i].i = i;
	std::cin >> s >> e;
	Polygon U, L;
	Pos u0 = e, u1 = e, u2 = s, u3 = s;
	Pos l0 = e, l1 = e, l2 = s, l3 = s;
	bool colu = 1, coll = 1;
	for (int i = 0; i < N; i++) {
		ll tq = cross(s, e, P[i]), fc = dot(s, e, P[i]);
		if (tq > 0 || (tq == 0 && fc > 0)) {
			U.push_back(P[i]);
			if (tq) colu = 0;
			if (dot(s, e, u0) > fc || (dot(s, e, u0) == fc && cross(s, e, u0) < tq)) u0 = P[i];
			if (dot(s, e, u1) > fc || (dot(s, e, u1) == fc && cross(s, e, u1) > tq)) u1 = P[i];
			if (dot(s, e, u2) < fc || (dot(s, e, u2) == fc && cross(s, e, u2) > tq)) u2 = P[i];
			if (dot(s, e, u3) < fc || (dot(s, e, u3) == fc && cross(s, e, u3) < tq)) u3 = P[i];
		}
		if (tq < 0 || (tq == 0 && fc < 0)) {
			L.push_back(P[i]);
			if (tq) coll = 0;
			if (dot(s, e, l0) > fc || (dot(s, e, l0) == fc && cross(s, e, l0) > tq)) l0 = P[i];
			if (dot(s, e, l1) > fc || (dot(s, e, l1) == fc && cross(s, e, l1) < tq)) l1 = P[i];
			if (dot(s, e, l2) < fc || (dot(s, e, l2) == fc && cross(s, e, l2) < tq)) l2 = P[i];
			if (dot(s, e, l3) < fc || (dot(s, e, l3) == fc && cross(s, e, l3) > tq)) l3 = P[i];
		}
	}
	if (U.empty() || L.empty()) {
		std::sort(P.begin(), P.end());
		std::cout << "YES\n"; 
		for (int i = 0; i < N - 1; i++) std::cout << P[i].i << " " << P[i + 1].i << "\n";
		return;
	}
	if (colu && coll) { std::cout << "NO\n"; return; }
	Pos u, l;
	bool ok = 0;
	int sz;
	if (cross(l1, u1, s) <= 0) { ok = 1; u = u1; l = l1; }
	else if (cross(l2, u2, e) >= 0) { ok = 1; u = u2; l = l2; }
	else if (u0 != u1 && l0 != l1 && cross(l0, u0, s) <= 0 && cross(l0, u0, u1) < 0 && cross(l0, u0, l1) < 0) { ok = 1; u = u0; l = l0; }
	else if (u3 != u2 && l3 != l2 && cross(l3, u3, e) >= 0 && cross(l3, u3, u2) > 0 && cross(l3, u3, l2) > 0) { ok = 1; u = u3; l = l3; }
	if (!ok) { std::cout << "NO\n"; return; }
	sz = U.size();
	std::cout << "YES\n" << u.i << " " << l.i << "\n";
	for (int i = 0; i < sz; i++) if (U[i] == u) { std::swap(U[i], u); break; }
	std::sort(U.begin() + 1, U.end(), [&](const Pos& p, const Pos& q) -> bool {
		ll tq = cross(U[0], p, q);
		return !tq ? p.Euc() < q.Euc() : tq > 0;
		});
	for (int i = 0; i < sz; i++) std::cout << U[i].i << " " << U[i + 1].i << "\n";
	sz = L.size();
	for (int i = 0; i < sz; i++) if (L[i] == l) { std::swap(L[i], l); break; }
	std::sort(L.begin() + 1, L.end(), [&](const Pos& p, const Pos& q) -> bool {
		ll tq = cross(L[0], p, q);
		return !tq ? p.Euc() < q.Euc() : tq > 0;
		});
	for (int i = 0; i < sz; i++) std::cout << L[i].i << " " << L[i + 1].i << "\n";
	return;
}
int main() { solve(); return 0; }