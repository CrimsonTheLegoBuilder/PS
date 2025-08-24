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
const ll MOD = 1e9 + 7;

int Q, N, M, K, a, b, c, d, T, W, H, L, R, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
int F[105];
int X;
//ld p, p0, X, Y, D, A, B, F;
//std::string S;
int S;
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
	std::cout.precision(1);
	std::cin >> N;
	while (N--) {
		std::cin >> M;
		K = M % 10;
		M /= 10;
		T += pow_mod(M, K);
	}
	std::cout << T << "\n";
 	return 0;
}