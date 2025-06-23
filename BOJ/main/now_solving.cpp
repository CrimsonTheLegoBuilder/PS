//refer to jiangly at QOJ
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <map>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
const ll INF = 1e17;
const int LEN = 1e5 + 10;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline ll sq(const ll& x) { return x * x; }

//freopen("../../../input_data/triathlon_tests/triath.20", "r", stdin);
//freopen("../../../input_data/triathlon_tests/triathlon_out.txt", "w", stdout);

#define TM 100'000ll

struct Pair { int i, j; };
int N, M, K, T, Q;
struct Sphere {
	int x, y, z, r;
	Sphere(int x_ = 0, int y_ = 0, int z_ = 0, int r_ = 0) : x(x_), y(y_), z(z_), r(r_) {}
	bool operator < (const Sphere& q) const { return r < q.r; }
	Sphere operator - (const Sphere& q) const { return Sphere(x - q.x, y - q.y, z - q.z, r + q.r); }
	Sphere& operator *= (const int& n) { x *= n, y *= n, z *= n, r *= n; return *this; }
	ll Euc() const { return (ll)x * x + (ll)y * y + (ll)z * z; }
	ll dist() const { return (ll)std::max(0., sqrt(Euc() - r) * .5); }
} S[LEN];
bool meet(const Sphere& p, const Sphere& q) {
	Sphere z = p - q;
	return z.Euc() - sq(z.r) <= 0;
}
bool meet(const int& i, const int& j) { return meet(S[i], S[j]); }
std::map<ll, Vint> MP;
//Vpii V;
std::vector<Pair> V;
ll hash(const int& x, const int& y, const int& z) {
	if (x < 0 || y < 0 || z < 0) return -1;
	if (x > TM || y > TM || z > TM) return -1;
	return x * TM * TM + y * TM + z;
}
void recur(const int& n) {
	if (n < 0) return;
	MP.clear();
	int d = S[n].r + 1 >> 1;
	int m = 0;
	while (S[m].r < d) m++;
	d <<= 2;
	for (int i = n; i >= 0; i--) {
		int x = S[i].x / d;
		int y = S[i].y / d;
		int z = S[i].z / d;
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				for (int dz = -1; dz <= 1; dz++) {
					ll hs = hash(x + dx, y + dy, z + dz);
					if (hs < 0) continue;
					auto it = MP.find(hs);
					if (it == MP.end()) continue;
					const Vint& I = it->second;
					for (const int& j : I) {
						if (meet(i, j)) V.emplace_back(i, j);
						if (V.size() >= K) return;
					}
				}
			}
		}
		if (i >= m) MP[hash(x, y, x)].push_back(i);
	}
	recur(m - 1);
	return;
}
int count(const int& m) {
	V.clear();
	for (int i = 0; i < N; i++) S[i].r += m - 1;
	recur(N - 1);
	for (int i = 0; i < N; i++) S[i].r -= m - 1;
	return V.size();
}
void query() {
	std::cin >> N >> K;
	for (int i = 0; i < N; i++)
		std::cin >> S[i].x >> S[i].y >> S[i].z >> S[i].r,
		S[i] *= 2;
	std::sort(S, S + N);
	int s = 0, e = 200'000;
	while (s < e) {
		int m = s + e >> 1;
		int cnt = count(m);
		if (cnt >= K) e = m;
		else s = m + 1;
	}
	if (s) {
		V.clear();
		for (int i = 0; i < N; i++) S[i].r += s - 1;
		recur(N - 1);
		for (int i = 0; i < N; i++) S[i].r -= s - 1;
	}
	Vll ret;
	for (const Pair& p : V) ret.push_back((S[p.i] - S[p.j]).dist());
	std::sort(ret.begin(), ret.end());
	int sz = ret.size();
	for (const ll& d : ret) std::cout << d << "\n";
	K -= sz;
	while (K--) std::cout << s << "\n";
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> Q; while (Q--) query();
	return;
}
int main() { solve(); return 0; }//boj18592
//refer to jiangly at QOJ

/*

1
8
2 2
-2 -2
2 -2
-2 2
3 0
-3 0
0 3
0 -3

1
10
-2 -4
5 2
7 7
10 -2
9 -5
-5 10
1 -5
4 -9
5 1
7 -7

1
4
4 5
-2 4
1 4
-5 -2

1
5
1 1
2 4
1 4
0 4
-1 1

3
4
0 3
1 3
3 1
3 0
4
-4 0
5 3
0 -4
-1 0
5
4 4
5 0
3 3
3 2
-4 2
5
1 1
2 4
1 4
0 4
-1 1

2
4
-1 -1
1 1
1 -1
-1 1
5
-1 -1
1 1
1 -1
-1 1
0 2


*/

#define BOJ
#ifdef BOJ
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
//int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068