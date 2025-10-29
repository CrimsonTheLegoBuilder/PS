#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef long long ll;
typedef double ld;
const ll INF = 1e17;

int N;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pos& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
} H[100000];
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
void init() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(8);
	std::cin >> N;
	for (int i = 0; i < N; i++) std::cin >> H[i].x >> H[i].y;
	return;
}
void rotating_calipers() {
	init();
	std::swap(H[0], *std::min_element(H, H + N));
	std::sort(H + 1, H + N, [&](const Pos& p, const Pos& q) {
		ll ret = cross(H[0], p, q);
		if (!ret) return (H[0] - p).Euc() < (H[0] - q).Euc();
		return ret > 0;
		});
	int S = -1;
	for (int i = 0; i < N; i++) {
		while (S >= 1 && cross(H[S - 1], H[S], H[i]) <= 0) S--;
		H[++S] = H[i];
	}
	N = S + 1;
	if (N == 2) { std::cout << (H[0] - H[1]).mag() * 2 << "\n"; return; }
	auto jaw = [&](const int& i, const int& f) -> ll {
		return (H[(i + 1) % N] - H[i]) / (H[(f + 1) % N] - H[f]);
		};
	auto w = [&](const int& i, const int& f) -> ld {
		const Pos& p0 = H[i], & p1 = H[(i + 1) % N], p2 = H[f];
		ll tq = cross(p0, p1, p2);
		return tq / (p1 - p0).mag();
		};
	ld ret = INF;
	for (int i = 0, j = 1; i < N; i++) {
		while (jaw(i, j) > 0) {
			j = (j + 1) % N;
		}
		ret = std::min(ret, w(i, j));
	}
	std::cout << ret << "\n";
	return;
}
int main() { rotating_calipers(); return 0; }