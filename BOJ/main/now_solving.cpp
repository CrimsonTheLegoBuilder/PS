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

int Q, N, M, Y, K, P, C, L, T; ll a, d, c, g, b, s, t, x, y;
int q1, q2, q3, q4, q0;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> T;
	while (T--) {
		std::vector<int> I(7);
		for (int& i : I) std::cin >> i;
		std::sort(I.begin(), I.end());
		N = 0; M = 0;
		for (const int& i : I) if (!(i & 1)) { M = i; break; }
		for (const int& i : I) if (!(i & 1)) N += i;
		std::cout << N << " " << M << "\n";
	}
	return 0;
}