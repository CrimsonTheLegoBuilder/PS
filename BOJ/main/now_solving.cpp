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
	std::cin >> A >> B >> C;
	F[0] = A;
	F[1] = B;
	F[2] = C;
	F[3] = A * B;
	F[4] = B * C;
	F[5] = C * A;
	F[6] = A * B * C;
	std::sort(F, F + 7);
	for (int i = 6; i > -1; i--) {
		if (F[i] & 1) { std::cout << F[i] << "\n"; return 0; }
	}
	std::cout << F[6] << "\n";
	return 0;
}