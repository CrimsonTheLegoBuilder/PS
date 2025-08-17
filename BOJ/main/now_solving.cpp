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

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, n[3], g,t, A, B, C, D, E, P, I;
int F[2001];
ld p, p0, X, Y;
std::string S;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> T;
	while (T--) {
		std::cin >> N;
		K = 0;
		for (int _ = 0; _ < N; _++) {
			std::cin >> M; K += M;
		}
		std::cout << (K % N ? "NO\n" : "YES\n");
	}
	return 0;
}