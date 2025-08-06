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

int Q, N, M, Y, K, P, C, L, T;
ll a, d, c, g, b, s, t, x, y, w1, h1, w2, h2;
int q1, q2, q3, q4, q0;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> w1 >> h1 >> w2 >> h2;
	T = 0;
	T += h1 + h2;
	T += std::max(w1, w2);
	std::cout << (T << 1) + 4 << "\n";
	return 0;
}