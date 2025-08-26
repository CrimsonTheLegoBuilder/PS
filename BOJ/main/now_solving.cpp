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
	std::cout.precision(3);
	std::cin >> N;
	while (N--) {
		ld a, b, c;
		std::cin >> a >> b >> c;
		ld a1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
		ld a2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
		if (a1 < a2) std::swap(a1, a2);
		std::cout << a1 << " " << a2 << "\n";
	}
	return 0;
}