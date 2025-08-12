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

int Q, N, M, K, a, b, T, W, H;
int q1, q2, q3, q4, q0, n[3], g, p ,t;
char B[1001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N >> B;
	for (int i = 0; i < N; i++) {
		if (B[i] == 'N') H++;
		if (B[i] == 'S') H--;
		if (B[i] == 'E') W++;
		if (B[i] == 'W') W--;

	}
	std::cout << std::abs(W) + std::abs(H) << "\n";
	return 0;
}