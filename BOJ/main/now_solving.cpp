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

int Q, N, M, K, a, b, c, T, W, H;
int q1, q2, q3, q4, q0, n[3], g, p ,t;
char B[1001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> Q;
	while (Q--) {
		int V[7], C = 3;
		memset(V, 0, sizeof V);
		while (C--) {
			std::cin >> N;
			V[N]++;
		}
		K = 0;
		for (int i = 1; i <= 6; i++) {
			if (V[i] == 3) {
				T = std::max(T, 10000 + i * 1000);
				K = 0;
				break;
			}
			if (V[i] == 2) {
				T = std::max(T, 1000 + i * 100);
				K = 0;
				break;
			}
			K = i;
		}
		if (K) T = std::max(T, K * 100);
	}
	std::cout << T << "\n";
	return 0;
}