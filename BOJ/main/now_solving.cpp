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
const ld PI = acos(-1);

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, S, I;
int F[8];
//ld p, p0, X, Y, D, A, B, F;
//std::string S;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N;
	F[1] = 1;
	while (N--) {
		std::cin >> a >> b;
		std::swap(F[a], F[b]);
	}
	for (int i = 1; i <= 3; i++) {
		if (F[i]) std::cout << i << "\n";
	}
	return 0;
}