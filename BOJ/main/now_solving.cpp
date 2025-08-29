#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <stdio.h>
typedef long long ll;
typedef long double ld;
//typedef double ld;
const ld PI = acos(-1);
const ld TOL = 1e-9;
const ll MOD = 1e9 + 7;
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
bool eq(const ld& x, const ld& y) { return !sign(x - y); }
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }

int Q, N, M, K, a, b, c, d, T, W, H, L, R, o, V;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
ll F[1000005];
int X;
//ld p, p0, X, Y, D, A, B, F;
std::string S, S0;
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
	std::cout.precision(9);
	int x0, y0, x, y, x1 = 0, y1 = 0;
	std::cin >> x0 >> y0 >> N;
	D = 1e9;
	while (N--) {
		std::cin >> x >> y;
		int d = std::abs(x - x0) + std::abs(y - y0);
		if (d < D) {
			D = d;
			x1 = x, y1 = y;
		}
	}
	std::cout << x1 << " " << y1 << "\n";
	return 0;
}