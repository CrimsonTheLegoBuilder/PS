#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
#include <deque>
#include <random>
#include <array>
#include <tuple>
#include <complex>
#include <numeric>
#include <set>
#include <chrono>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
const ld INF = 1e18;
const ld TOL = 1e-9;
const ld PI = acos(-1);
const int LEN = 1e4;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ll sq(const ll& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo = 0, const ld& hi = 1) { return std::min(hi, std::max(lo, x)); }
inline ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
inline ll gcd(ll x, ll y, ll z) {
	x = std::abs(x); y = std::abs(y); z = std::abs(z);
	ll w = gcd(x, y);
	return gcd(w, z);
}

std::mt19937 gen(std::chrono::steady_clock::now().time_since_epoch().count());

bool DEBUG = 0;

#define LOCAL_TEST

const int N_ = 8000;
const int K_ = 140;
const int DX = 14;
const int DY = 10;

int N, K;
struct Pii {
	int x, y; int i;
	Pii(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pii& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pii& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pii& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pii& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pii operator + (const Pii& p) const { return { x + p.x, y + p.y }; }
	Pii operator - (const Pii& p) const { return { x - p.x, y - p.y }; }
	Pii operator * (const int& n) const { return { x * n, y * n }; }
	Pii operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pii& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pii& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	Pii& operator += (const Pii& p) { x += p.x; y += p.y; return *this; }
	Pii& operator -= (const Pii& p) { x -= p.x; y -= p.y; return *this; }
	Pii& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pii& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pii operator - () const { return { -x, -y }; }
	Pii operator ~ () const { return { -y, x }; }
	Pii operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ll Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pii& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pii& p) { os << p.x << " " << p.y; return os; }
};
const Pii Oii = { 0, 0 };
const Pii INF_PT = { (int)1e9, (int)1e9 };
typedef std::vector<Pii> Vpii;
bool cmpx(const Pii& p, const Pii& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pii& p, const Pii& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpi(const Pii& p, const Pii& q) { return p.i < q.i; }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return sign(cross(d1, d2, d3, d4)); }
struct Order {
	ll o;
	int i;
	bool operator < (const Order& q) const { return o < q.o; }
};
ll hilbert_order(const Pii& p, int pow2) {
	int x = p.x, y = p.y;
	ll d = 0;
	for (int s = pow2 >> 1; s > 0; s >>= 1) {
		int rx = (x & s) > 0;
		int ry = (y & s) > 0;
		d += (ll)s * s * ((3ll * rx) ^ ry);
		if (ry == 0) {
			if (rx == 1) {
				x = pow2 - 1 - x;
				y = pow2 - 1 - y;
			}
			int temp = x; x = y; y = temp;
		}
	}
	return d;
}
Vpii clst[DX][DY];
int clst_cnt[K_], clst_cnt_2d[DX][DY];
int fst_belt_cnt[DX];
int grp[LEN];
Vint dt[LEN];
bool F[K_];
struct ClusterMeta { ld d; ll sx, sy; } meta[DX][DY];
void update_meta(int x, int y) {
	const Vpii& path = clst[x][y];
	int sz = path.size();
	if (sz == 0) { meta[x][y] = { 0.0, 0, 0 }; return; }
	ld d = 0;
	ll sx = 0, sy = 0;
	for (int i = 0; i < sz; i++) {
		const Pii& p0 = path[i], & p1 = path[(i + 1) % sz];
		d += (p0 - p1).mag();
		sx += p0.x; sy += p0.y;
	}
	meta[x][y] = { d, sx, sy };
	return;
}
void init_all_meta() { for (int x = 0; x < DX; x++) for (int y = 0; y < DY; y++) update_meta(x, y); }
Pii centroid(int x, int y) {
	ll sx = 0, sy = 0;
	int sz = clst[x][y].size();
	if (sz == 0) return { -1, -1 };
	for (const Pii& p : clst[x][y]) { sx += p.x; sy += p.y; }
	return { (int)(sx / sz), (int)(sy / sz) };
}
Pii centroid_fast(int x, int y) {
	int sz = clst[x][y].size();
	if (sz == 0) return { -1, -1 };
	return { (int)(meta[x][y].sx / sz), (int)(meta[x][y].sy / sz) };
}
void first_clustering(Vpii& P) {
	for (int i = 0; i < K; i++) {
		if (!((i + 1) % 7)) clst_cnt[i] = clst_cnt_2d[i / DY][i % DY] = 58;
		else clst_cnt[i] = clst_cnt_2d[i / DY][i % DY] = 57;
	}
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			fst_belt_cnt[x] += clst_cnt_2d[x][y];
		}
	}
	int prv_x = 0;
	int g = 0;
	memset(grp, -1, sizeof grp);
	std::sort(P.begin(), P.end(), cmpx);
	for (int x = 0; x < DX; x++) {
		Vpii C;
		for (int i = 0; i < fst_belt_cnt[x]; i++) C.push_back(P[prv_x + i]);
		std::sort(C.begin(), C.end(), cmpy);
		int prv_y = 0;
		for (int y = 0; y < DY; y++) {
			for (int i = 0; i < clst_cnt_2d[x][y]; i++) {
				grp[C[prv_y + i].i] = g;
				clst[x][y].push_back(C[prv_y + i]);
			}
			g++;
			prv_y += clst_cnt_2d[x][y];
		}
		prv_x += fst_belt_cnt[x];
	}
	std::sort(P.begin(), P.end(), cmpi);
	return;
}
void init_hilbert_paths(Vpii& P) {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			if (clst[x][y].empty()) continue;
			std::vector<Order> V;
			int min_x = 1e9, min_y = 1e9;
			int max_x = -1e9, max_y = -1e9;
			int sz = clst[x][y].size();
			for (int i = 0; i < sz; i++) {
				min_x = std::min(min_x, clst[x][y][i].x);
				min_y = std::min(min_y, clst[x][y][i].y);
				max_x = std::max(max_x, clst[x][y][i].x);
				max_y = std::max(max_y, clst[x][y][i].y);
			}
			ll pow2 = 1;
			int tx = max_x - min_x;
			int ty = max_y - min_y;
			int t = std::max(tx, ty);
			while (pow2 <= t) pow2 <<= 1;
			for (int i = 0; i < clst_cnt_2d[x][y]; i++) {
				Pii p = clst[x][y][i] - Pii(min_x, min_y);
				V.push_back({ hilbert_order(p, pow2), clst[x][y][i].i });
			}
			std::sort(V.begin(), V.end());
			Vpii hlbrt;
			hlbrt.reserve(sz);
			for (int i = 0; i < sz; i++) hlbrt.push_back(P[V[i].i]);
			clst[x][y] = hlbrt; // µ¤¾î¾²±â
		}
	}
	init_all_meta();
	return;
}
void optimize_2opt_single(int x, int y, int thr = 100000) {
	Vpii& path = clst[x][y];
	int sz = path.size();
	if (sz < 3) return;
	bool imp = 1;
	while (imp && thr--) {
		imp = 0;
		for (int i = 0; i < sz - 1; i++) {
			for (int j = i + 2; j < sz; j++) {
				if (i == 0 && j == sz - 1) continue;
				Pii p1 = path[i];
				Pii p2 = path[i + 1];
				Pii p3 = path[j];
				Pii p4 = path[(j + 1) % sz];
				ld d12 = (p1 - p2).mag();
				ld d34 = (p3 - p4).mag();
				ld d13 = (p1 - p3).mag();
				ld d24 = (p2 - p4).mag();
				if (d13 + d24 < d12 + d34 - TOL) {
					std::reverse(path.begin() + i + 1, path.begin() + j + 1);
					imp = 1;
				}
			}
		}
	}
	update_meta(x, y);
	return;
}
void optimize_2opt(int thr = 100000) {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			optimize_2opt_single(x, y, thr);
		}
	}
	return;
}
void apply_double_bridge(int x, int y) {
	Vpii& path = clst[x][y];
	int sz = path.size();
	if (sz < 8) return;

	std::vector<int> cuts;
	std::uniform_int_distribution<int> dist_idx(1, sz - 1);

	while (cuts.size() < 3) {
		int k = dist_idx(gen);
		bool dup = 0;
		for (int c : cuts) if (c == k) dup = 1;
		if (!dup) cuts.push_back(k);
	}
	std::sort(cuts.begin(), cuts.end());
	int A = cuts[0];
	int B = cuts[1];
	int C = cuts[2];
	//1-4-3-2
	Vpii new_path;
	new_path.reserve(sz);
	for (int i = 0; i < A; i++) new_path.push_back(path[i]);
	for (int i = C; i < sz; i++) new_path.push_back(path[i]);
	for (int i = B; i < C; i++) new_path.push_back(path[i]);
	for (int i = A; i < B; i++) new_path.push_back(path[i]);
	clst[x][y] = new_path;
	return;
}
void optimize_2opt_kick(int thr = 100000) {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			optimize_2opt_single(x, y, thr);
			int kick_iter = 100;
			ld best_d = meta[x][y].d;
			Vpii best_path = clst[x][y];
			while (kick_iter--) {
				Vpii backup = clst[x][y];
				apply_double_bridge(x, y);
				optimize_2opt_single(x, y, 100000);
				ld cur_d = meta[x][y].d;
				if (cur_d < best_d - 1e-5) {
					best_d = cur_d;
					best_path = clst[x][y];
				}
				else {
					clst[x][y] = backup;
					update_meta(x, y);
				}
			}
			clst[x][y] = best_path;
		}
	}
	//init_all_meta();
	return;
}
void reset_data() {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			clst[x][y].clear();
		}
	}
	memset(clst_cnt, 0, sizeof(clst_cnt));
	memset(clst_cnt_2d, 0, sizeof(clst_cnt_2d));
	memset(fst_belt_cnt, 0, sizeof(fst_belt_cnt));
	memset(grp, 0, sizeof(grp));
	return;
}
ld run_solver() {
	std::cin >> N >> K;
	Vpii P(N);
	for (Pii& p : P) std::cin >> p;
	for (int i = 0; i < N; i++) P[i].i = i;
	first_clustering(P);
	init_hilbert_paths(P);
	optimize_2opt();
	optimize_2opt_kick();
	ld max_dist = 0;
	int min_idx = 1e9, max_idx = -1;
	std::set<int> S;
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			if (clst[x][y].empty()) continue;
#ifndef LOCAL_TEST
			std::cout << clst[x][y].size() << "\n";
			for (const Pii& p : clst[x][y]) std::cout << p.i + 1 << " ";
			std::cout << "\n";
#else
			ld cur_dist = 0;
			int sz = clst[x][y].size();
			const Vpii& H = clst[x][y];
			for (int i = 0, u, v; i < sz; i++) {
				u = H[i].i; v = H[(i + 1) % sz].i;
				min_idx = std::min(min_idx, u);
				max_idx = std::max(max_idx, u);
				S.insert(u);
				cur_dist += (H[i] - H[(i + 1) % sz]).mag();
			}
			max_dist = std::max(max_dist, cur_dist);
#endif
		}
	}
#ifdef LOCAL_TEST
	std::cout << "DEBUG::\n";
	std::cout << "min:: " << min_idx << "\n";
	std::cout << "max:: " << max_idx << "\n";
	std::cout << "sz :: " << S.size() << "\n";
	std::cout << "DEBUG::\n";
#endif
	return max_dist;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(10);
#ifdef LOCAL_TEST
	ld total_score = 0;
	ld max_score = 0;
	int tc = 50;
	std::cout << "========= [LOCAL TEST START] =========\n";
	for (int i = 1; i <= tc; i++) {
		std::string filename = (i < 10 ? "0" : "") + std::to_string(i) + ".in";
		std::string path = "../../tests/814_3/" + filename;
		if (freopen(path.c_str(), "r", stdin) == NULL) { std::cout << "File Not Found: " << path << "\n"; continue; }
		reset_data();
		ld score = run_solver();
		std::cout << "Test " << filename << " : " << score << "\n";
		total_score += score;
		max_score = std::max(max_score, score);
	}
	std::cout << "======================================\n";
	std::cout << "2-opt:\n";
	std::cout << "Total Score: " << total_score << "\n";
	std::cout << "Average Score: " << total_score / 50.0 << "\n";
	std::cout << "Worst Case   : " << max_score << "\n";
	std::cout << "======================================\n";
#else
	reset_data();
	run_solver();
#endif
	return;
}
int main() { solve(); return 0; }
