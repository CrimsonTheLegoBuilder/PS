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
int F[200005];
ll TT[305];
ll TTT[305];
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
	std::cout.precision(14);
	std::cin >> N;
	T = 0;
	for (int i = 0; i < N; i++) {
		std::cin >> K;
		K = std::max(K - (N - i), 0);
		T = std::max(T, K);
	}
	std::cout << T << "\n";
	return 0;
}