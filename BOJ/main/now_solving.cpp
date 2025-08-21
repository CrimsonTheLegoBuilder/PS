#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
typedef long long ll;
//typedef long double ld;
typedef double ld;
const ld PI = acos(-1);
const ll MOD = 998244353;

int Q, N, M, K, a, b, c, d, T, W, H, L, R, s, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
int F[5];
int X;
//ld p, p0, X, Y, D, A, B, F;
std::string S[101];
ll pow_mod(ll N, int K) {
	if (!K) return 1;
	ll h = pow_mod(N, K / 2);
	ll ret = (h * h) % MOD;
	if (K & 1) ret = (ret * N) % MOD;
	return ret % MOD;
}
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	while (1) {
		std::cin >> A >> B >> C >> D;
		if (!A && !B && !C && !D) break;
		if (!A) std::cout << D / (B * C) << " " << B << " " << C << " " << D << "\n";
		if (!B) std::cout << A << " " << D / (A * C) << " " << C << " " << D << "\n";
		if (!C) std::cout << A << " " << B << " " << D / (A * B) << " " << D << "\n";
		if (!D) std::cout << A << " " << B << " " << C << " " << A * B * C << "\n";
	}
	return 0;
}