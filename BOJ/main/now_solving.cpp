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
	std::cout.precision(12);
	std::cin >> Q;
	int I[7] = { 2, 3, 4, 5, 6, 7, 2 };
	T = 0;
	for (int i = 0; i < 7; i++) {
		T += (Q % 10) * I[i];
		Q /= 10;
	}
	//N = T % 11;
	if (N == 0) std::cout << "J\n";
	if (N == 1) std::cout << "A\n";
	if (N == 2) std::cout << "B\n";
	if (N == 3) std::cout << "C\n";
	if (N == 4) std::cout << "D\n";
	if (N == 5) std::cout << "E\n";
	if (N == 6) std::cout << "F\n";
	if (N == 7) std::cout << "G\n";
	if (N == 8) std::cout << "H\n";
	if (N == 9) std::cout << "I\n";
	if (N == 10) std::cout << "Z\n";
	return 0;
}