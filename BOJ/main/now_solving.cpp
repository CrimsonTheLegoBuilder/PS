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

int Q, N, M, K, a, b, c, T, W, H, s, o;
int q1, q2, q3, q4, q0, n[3], g, p ,t, A, I;
char B[2001][2001];
int ret[2001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N >> M;
	ret[1] = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			std::cin >> B[i][j];
			if (B[i][j] == 'X') ret[1]++;
		}
	}
	for (int k = 1; k < std::min(N, M); k++) {
		//std::cout << "k:: " << k << "\n";
		for (int i = 0; i < N - k; i++) {
			for (int j = 0; j < M - k; j++) {
				bool f = 1;
				for (int l = 0; l <= k; l++) {
					if (B[i + l][j + l] != 'X') f = 0;
				}
				if (f) {
					int f0 = 1;
					for (int n = i; n <= i + k; n++) {
						for (int m = j; m <= j + k; m++) {
							if (n - i == m - j) continue;
							if (B[n][m] != '.') {
								f0 = 0;
								break;
							}
						}
						if (!f0) break;
					}
					if (f0) ret[k + 1]++;
				}
			}
		}
	}
	for (int i = 1; i <= std::min(N, M); i++) {
		std::cout << ret[i] << "\n";
	}
	return 0;
}