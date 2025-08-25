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

int Q, N, M, K, a, b, c, d, T, W, H, L, R, o, V;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
ll F[1000005];
int X;
//ld p, p0, X, Y, D, A, B, F;
std::string S;
//int S;
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
	std::cin >> T;
	ll J = 0, j = 0, m = 0;
	while (T--) {
		std::cin >> S >> N;
		if (S == "J") j++, J += N;
		if (S == "M") m++, M += N;
	}
	if (J * m == M * j) std::cout << "V\n";
	else if (J * m > M * j) std::cout << "J\n";
	//else if (J * m < M * j) std::cout << "M\n";
	return 0;
}