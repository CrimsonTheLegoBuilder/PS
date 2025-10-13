#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
typedef long long ll;

int N;
ll A, B;
struct Pos { ll x, y; } P[1 << 7];
Pos dir[4] = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
int inner(const Pos& p, const ll& x1, const ll& y1, const ll& x2, const ll& y2) {
	return (x1 <= p.x && y1 <= p.y && p.x <= x2 && p.y <= y2) ? 1 : 0;
}
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);
	std::cin >> N >> A >> B;
	int ret = 0;
	for (int i = 0; i < N; i++) std::cin >> P[i].x >> P[i].y;
	for (int i = 0; i < N; i++) {
		for (int d = 0; d < 4; d++) {
			int c = 0;
			ll x0 = P[i].x + A * dir[d].x;
			ll y0 = P[i].y + B * dir[d].y;
			ll x1 = std::min(P[i].x, x0);
			ll y1 = std::min(P[i].y, y0);
			ll x2 = std::max(P[i].x, x0);
			ll y2 = std::max(P[i].y, y0);
			for (int j = 0; j < N; j++) {
				c += inner(P[j], x1, y1, x2, y2);
			}
			ret = std::max(ret, c);
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (std::abs(P[i].x - P[j].x) > A) continue;
			if (std::abs(P[i].y - P[j].y) > B) continue;
			ll x00 = std::min(P[i].x, P[j].x);
			ll y00 = std::min(P[i].y, P[j].y);
			ll x11 = std::max(P[i].x, P[j].x);
			ll y11 = std::max(P[i].y, P[j].y);
			ll x[2] = { x00, x11 };
			ll y[2] = { y00, y11 };
			for (int u = 0; u < 2; u++) {
				for (int v = 0; v < 2; v++) {
					for (int d = 0; d < 4; d++) {
						int c = 0;
						ll x0 = x[u] + A * dir[d].x;
						ll y0 = y[v] + B * dir[d].y;
						ll x1 = std::min(x[u], x0);
						ll y1 = std::min(y[v], y0);
						ll x2 = std::max(x[u], x0);
						ll y2 = std::max(y[v], y0);
						for (int k = 0; k < N; k++) {
							c += inner(P[k], x1, y1, x2, y2);
						}
						ret = std::max(ret, c);
					}
				}
			}
		}
	}
	std::cout << ret << "\n";
	return 0;
}