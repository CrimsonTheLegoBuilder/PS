#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef long long ll;
typedef double ld;
const ll INF = 1e18;
const int LEN = 3e4 + 10;

int N, T;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	ll operator / (const Pos& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
} P[LEN], V[LEN], H[LEN], C[LEN];
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll rotating_calipers(const int& t) {
	//for (int i = 0; i < N; i++) H[i] = P[i] + V[i] * t;
	//std::swap(H[0], *std::min_element(H, H + N));
	//std::sort(H + 1, H + N, [&](const Pos& p, const Pos& q) {
	//	ll ret = cross(H[0], p, q);
	//	return ret == 0 ? (H[0] - p).Euc() < (H[0] - q).Euc() : ret > 0;
	//	});
	//int S = -1;
	//for (int i = 0; i < N; i++) {
	//	while (S >= 1 && cross(H[S - 1], H[S], H[i]) <= 0) S--;
	//	H[++S] = H[i];
	//}
	for (int i = 0; i < N; i++) C[i] = P[i] + V[i] * t;
	std::sort(C, C + N);
	int S = -1;
	for (int i = 0; i < N; i++) {
		while (S > 0 && cross(H[S - 1], H[S], C[i]) <= 0) S--;
		H[++S] = C[i];
	} S--;
	int S1 = S + 1;
	for (int i = N - 1; i >= 0; i--) {
		while (S > S1 && cross(H[S - 1], H[S], C[i]) <= 0) S--;
		H[++S] = C[i];
	} S--;
	int sz = S + 1;
	if (sz == 1) return 0;
	if (sz == 2) return (H[0] - H[1]).Euc();
	ll d = 0;
	for (int i = 0, j = 1; i < sz; i++) {
		while (cross(H[i], H[(i + 1) % sz], H[j], H[(j + 1) % sz]) >= 0) {
			d = std::max(d, (H[i] - H[j]).Euc());
			j = (j + 1) % sz;
		}
		d = std::max(d, (H[i] - H[j]).Euc());
	}
	return d;
}
ll ternary_search() {
	ll s = 0, e = T, d1, d2;
	while (s + 2 < e) {
		ll t1 = (s + s + e) / 3;
		ll t2 = (s + e + e) / 3;
		d1 = rotating_calipers(t1);
		d2 = rotating_calipers(t2);
		if (d1 <= d2) e = t2;
		else s = t1;
	}
	ll d = INF;
	for (int t = s; t <= e; t++) {
		ll l = rotating_calipers(t);
		if (d > l) d = l, T = t;
	}
	return d;
}
ll bi_search() {
	ll s = 0, e = T, d1, d2;
	while (s + 1 < e) {
		ll t = s + e >> 1;
		d1 = rotating_calipers(t);
		d2 = rotating_calipers(t + 1);
		if (d1 <= d2) e = t;
		else s = t;
	}
	ll d = INF;
	for (int t = s; t <= e; t++) {
		ll l = rotating_calipers(t);
		if (d > l) d = l, T = t;
	}
	return d;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N >> T;
	for (int i = 0; i < N; i++) std::cin >> P[i].x >> P[i].y >> V[i].x >> V[i].y;
	//ll D = ternary_search();
	ll D = bi_search();
	std::cout << T << "\n" << D << "\n";
	return;
}
int main() { solve(); return 0; }//boj13310