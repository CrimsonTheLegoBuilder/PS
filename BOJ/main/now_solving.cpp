#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <unordered_set>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::unordered_set<int> Sint;
const int LEN = 1e5;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }

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
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	int i = 0;
	std::cin >> N; Polygon P(N); for (Pos& p : P) std::cin >> p, p.i = i++;
	Sint F, B;
	Polygon S;
	for (int i = 0; i < N; i++) {
		while (S.size() > 1 && ccw(S[S.size() - 2], S[S.size() - 1], P[i]) < 0)
			S.pop_back();
		S.push_back(P[i]);
	}
	for (const Pos& p : S) F.insert(p.i);
	S.clear();
	for (int i = N - 1; i >= 0; i--) {
		while (S.size() > 1 && ccw(S[S.size() - 2], S[S.size() - 1], P[i]) > 0)
			S.pop_back();
		S.push_back(P[i]);
	}
	for (const Pos& p : S) B.insert(p.i);
	ld R = 0;
	for (int i = 0; i < N; i++) {
		if (F.find(i) != F.end() && B.find(i) != B.end()) {

		}
	}
	return;
}
int main() { solve(); return 0; }