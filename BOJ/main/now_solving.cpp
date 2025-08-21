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
int q1, q2, q3, q4, q0, A[100], B[100], C, D, E, P, I;
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
	std::cin >> N;
	if (N > 198) { std::cout << "0\n"; return 0; }
	for (int i = 0; i <= std::min(99, N); i++) {
		K = N - i;
		if (K < 100) T++;
	}
	std::cout << T << "\n";
	return 0;
}