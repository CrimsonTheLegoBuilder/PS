#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
//typedef double ld;
const ld PI = acos(-1);
const ld TOL = 1e-9;
const ll MOD = 1e9 + 7;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
bool eq(const ld& x, const ld& y) { return !sign(x - y); }

int Q, N, M, K, a, b, c, d, T, W, H, L, R, o, V;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
ll F[1000005];
int X;
//ld p, p0, X, Y, D, A, B, F;
//std::string S, S0;
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
	std::cout.precision(2);
	while (1) {
		std::cin >> N;
		if (!N) break;
		if (N <= 5) { std::cout << N << " DEFICIENT\n"; continue; }
		M = sqrt(N) + 1;
		T = 1;
		//std::cout << "M:: " << M << "\n";
		for (int i = 2; i < M; i++) {
			if (!(N % i) && i != (N / i)) T += i, T += N / i;
			else if (!(N % i)) T += i;
		}
		//std::cout << "N:: " << N << " T:: " << T << "\n";
		if (N == T) std::cout << "PERFECT\n";
		else if (N < T) std::cout << "ABUNDANT\n";
		else std::cout << "DEFICIENT\n";
	}
	return 0;
}