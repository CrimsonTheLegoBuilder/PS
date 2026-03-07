#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cassert>
const int LEN = 10000;
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }

bool B[LEN][LEN];
int N, sx, ex, sy, ey;
int SV[15], SH[15];
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
};
Pos dir[6] = { Pos(1, 0), Pos(1, 1), Pos(-1, 1), Pos(-1, 0), Pos(-1, -1), Pos(1, -1) };
void recur(Pos ps, Pos pe, int d, int s) {
	if (s == 1) {
		sx = std::min({ sx, ps.x, pe.x });
		ex = std::max({ ex, ps.x, pe.x });
		sy = std::min({ sy, ps.y, pe.y });
		ey = std::max({ ey, ps.y, pe.y });
		if (d == 0 || d == 3) {
			assert(ps.y == pe.y);
			int sx_ = std::min(ps.x, pe.x);
			int ex_ = std::max(ps.x, pe.x);
			int y = ps.y;
			for (int i = sx_; i <= ex_; i++) B[y][i] = 1;
		}
		else {
			int fx = std::abs(ps.x - pe.x);
			int fy = std::abs(ps.y - pe.y);
			assert(fx == fy);
			int dx = sign(pe.x - ps.x);
			int dy = sign(pe.y - ps.y);
			int x = ps.x, y = ps.y;
			while (fx--) {
				B[y][x] = 1;
				x += dx;
				y += dy;
			}
		}
		return;
	}

	Pos pv = dir[d];
	int itx = pv.y ? SV[s - 1] - 1 : SH[s - 1] - 1;
	int ity = pv.y ? SV[s - 1] - 1 : 0;
	pv.x *= itx; pv.y *= ity;
	Pos p0 = ps, p1 = ps + pv, p2 = ps + pv * 2, p3 = ps + pv * 3;
	assert(p3 == pe);

	recur(p0, p1, d, s - 1); recur(p2, p3, d, s - 1);

	pv = dir[(d + 5) % 6];
	itx = pv.y ? SV[s - 1] - 1 : SH[s - 1] - 1;
	ity = pv.y ? SV[s - 1] - 1 : 0;
	pv.x *= itx; pv.y *= ity;
	Pos pm = p1 + pv;

	recur(p1, pm, (d + 5) % 6, s - 1); recur(pm, p2, (d + 1) % 6, s - 1);

	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	//freopen("../tests/koch/ans/04.out", "w", stdout);
	std::cin >> N;

	if (N == 1) {
		std::cout << "  *\n";
		std::cout << " * *\n";
		std::cout << "*****\n";
		return;
	}

	SV[1] = 3; SH[1] = 5;
	for (int i = 2; i <= 11; i++) {
		SV[i] = SV[i - 1] * 3 - 2;
		SH[i] = SH[i - 1] * 3 - 2;
	}

	sx = LEN >> 1;
	sy = LEN >> 1;

	Pos p0, p1, p2;
	p0 = Pos(sx, sy) - Pos(SH[N] / 2, SV[N - 1] / 2);
	p1 = p0 + Pos(SH[N] - 1);
	p2 = p0 + Pos(SV[N] - 1, SV[N] - 1);

	recur(p0, p1, 0, N); recur(p1, p2, 2, N); recur(p2, p0, 4, N);

	for (int y = sy; y <= ey; y++) {
		int i = ex;
		for (i; i > 0; i--) if (B[y][i]) break;
		for (int x = sx; x <= i; x++) {
			std::cout << (B[y][x] ? "*" : " ");
		}
		std::cout << "\n";
	}
	//std::cout << tt << "\n";
	return;
}
int main() { solve(); return 0; }