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

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, n[3], g, p ,t, A, B, C, D, P, I;
int F[2001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> W >> H >> N >> A >> B;
	if (W < A || H < B) { std::cout << "-1\n"; return 0; }
	T = (W / A) * (H / B);
	std::cout << (N / T) + (!(N % T) ? 0 : 1) << "\n";
	return 0;
}