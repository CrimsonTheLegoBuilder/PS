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
	if (N == 1 || N == 2) std::cout << 4 << "\n";
	else if (N & 1) std::cout << N + 1;
	else std::cout << N << "\n";
	return 0;
}