#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
const int LEN = 1e5 + 10;

int N, T;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pos& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
} P[LEN], V[LEN], H[LEN];
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
void init() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(6);
	return;
}
ll rotating_calipers(const int& t) {
	init();
	std::swap(H[0], *std::min_element(H, H + N));
	std::sort(H + 1, H + N, [&](const Pos& p, const Pos& q) {
		ll ret = cross(H[0], p, q);
		return ret == 0 ? (H[0] - p).Euc() < (H[0] - q).Euc() : ret > 0;
		});
	int S = -1;
	for (int i = 0; i < N; i++) {
		while (S >= 1 && cross(H[S - 1], H[S], H[i]) <= 0) S--;
		H[++S] = H[i];
	}
	N = S + 1;
	if (N == 2) {  }
	return;
}
ll ternary_search() {
	ll s = 0, e = T, d1, d2;
	while (s + 2 < e) {
		ll t1 = (s + s + s + s + e + e + e) / 7;
		ll t2 = (s + s + s + e + e + e + e) / 7;
		d1 = rotating_calipers(t1);
		d2 = rotating_calipers(t2);
		if (d1 <= d2) e = t2;
		else s = t1;
	}
	ll d = INF;
	for (int t = s; t <= e; t++) d = std::min(d, rotating_calipers(t));
	return d;
}
#define BOJ
#ifdef BOJ
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N >> T;
	for (int i = 0; i < N; i++) std::cin >> P[i].x >> P[i].y >> V[i].x >> V[i].y;
	std::cout << ternary_search() << "\n";
	return;
}
#else
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	int Q = 130;
	for (int q = 0; q <= Q; q++) {
		std::string Din = "../../tests/aed/input/input";
		std::string Dout = "../../tests/aed/output/output";
		std::string Qin = Din + std::to_string(q) + ".txt";
		std::string Qout = Dout + std::to_string(q) + ".txt";
		freopen(Qin.c_str(), "r", stdin);
		//freopen("../../tests/aed/input/ret.txt", "w", stdout);
		Pos u;
		std::cin >> N; Polygon H(N); std::cin >> u;
		for (Pos& p : H) std::cin >> p; norm(H);
		Vint I, J;
		for (int i = 0, i1; i < N; i++) {
			i1 = (i + 1) % N;
			Pos p0 = H[i], p1 = H[i1];
			bool vis;
			if (!ccw(u, p0, p1)) {
				if (dot(p0, p1, u) < 0) std::swap(p0, p1);
				vis = 1;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					if (intersect(q0, q1, p1, u)) { vis = 0; break; }
				}
			}
			else {
				if (ccw(u, p0, p1) < 0) std::swap(p0, p1);
				Seg s1 = Seg(p0, p1);
				Vrange V;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					Seg s2 = Seg(q0, q1);
					Range r = range(s1, s2, u);
					if (r.s.den != -1) V.push_back(r);
				}
				vis = 0;
				std::sort(V.begin(), V.end());
				V.push_back(Range(Frac(1), Frac(1)));
				int sz = V.size();
				Frac hi = Frac(0);
				for (const Range& r : V) {
					if (hi < r.s) { vis = 1; break; }
					else hi = std::max(hi, r.e);
				}
			}
			if (vis) I.push_back(i + 1);
		}
		//std::cout << I.size() << "\n";
		//for (int i : I) std::cout << i << " ";
		freopen(Qout.c_str(), "r", stdin);
		std::cin >> M;
		J.resize(M); for (int& j : J) std::cin >> j;
		if (I.size() != M) { std::cout << "fuck::\n"; continue; }
		bool f = 1;
		for (int j = 0; j < M; j++) {
			if (I[j] != J[j]) { f = 0; std::cout << "fuck::\n"; break; }
		}
		if (f) std::cout << "good::\n";
		else std::cout << "fuck::\n";
	}
	return;
}
#endif
int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068
