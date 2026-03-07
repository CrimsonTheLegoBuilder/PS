#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <set>
//#define TIME
#ifdef TIME
#include <iomanip>
#include <chrono>
#endif
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 2e3 + 5;

//#define ROTATE

int W, H, N, M, Q;
int XL[LEN], XR[LEN];
int YL[LEN], YH[LEN];
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	Pos operator ~ () const { return { -y, x }; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} D0, D1;
const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon P[LEN];
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
};
bool check(const Polygon& p, const Polygon& b) {
	assert(p.size() == 2);
	if (p[0].x >= b[1].x ||
		p[1].x <= b[0].x ||
		p[0].y >= b[1].y ||
		p[1].y <= b[0].y)
		return 1;
	return 0;
}
bool check(const Polygon& b) {
	assert(b.size() == 2);
	for (int i = 1; i <= N; i++) {
		if (!check(P[i], b)) return 0;
	}
	return 1;
}
ll area(const Polygon& b) {
	assert(b.size() == 2);
	return ((ll)b[1].x - b[0].x) * ((ll)b[1].y - b[0].y);
}
void solve() {
#ifdef TIME
	auto start = std::chrono::high_resolution_clock::now();
#endif
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	//freopen("../tests/awesome/15.in", "r", stdin);
	//freopen("../tests/awesome/15_naive.out", "w", stdout);
	std::cin >> W >> H >> N;
	D0 = Pos(0, 0); D1 = Pos(W, H);
	P[0] = { D0, D1 };
	XL[0] = 0; XR[0] = W;
	YL[0] = 0; YH[0] = H;
	int xl, yd, xr, yh;
	for (int i = 1; i <= N; i++) {
		std::cin >> xl >> yd >> xr >> yh;
		Pos a, b;
		a = Pos(xl, yd); b = Pos(xr, yh);
		P[i] = { a, b };
		XL[i] = xr; XR[i] = xl;
		YL[i] = yh; YH[i] = yd;
	}
	ll ret = 0;
	for (int u = 0; u <= N; u++) {
		yh = YH[u];
		for (int r = 0; r <= N; r++) {
			xr = XR[r];
			for (int l = 0; l <= N; l++) {
				xl = XL[l];
				if (xr <= xl) continue;
				for (int d = 0; d <= N; d++) {
					yd = YL[d];
					if (yh <= yd) continue;
					Polygon B = { Pos(xl, yd), Pos(xr, yh) };
					if (check(B)) ret = std::max(ret, area(B));
				}
			}
		}
	}
	
	std::cout << ret << "\n";
#ifdef TIME
	std::cout << "\n\n\n";
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Start time: " << std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count() << " us\n";
	std::cout << "End time: " << std::chrono::duration_cast<std::chrono::microseconds>(end.time_since_epoch()).count() << " us\n";
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Execution time: " << duration.count() << " microseconds\n";
#endif
	return;
}
int main() { solve(); return 0; }//kitpc? 1X? CSS IS AWESOME