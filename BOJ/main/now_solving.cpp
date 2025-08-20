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

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
int F[8];
ll TT[305];
ll TTT[305];
//ld p, p0, X, Y, D, A, B, F;
std::string S[101];
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(14);
	std::cin >> N;
	while (N--) {
		std::cin >> M;
		std::cout << "Pairs for ";
		std::cout << M << ": ";
		K = M - 1;
		for (int i = 1; i <= K / 2; i++) {
			std::cout << i << " " << M - i;
			if (i <= K / 2 - 1) std::cout << ", ";
			else std::cout << "\n";
		}
	}
	return 0;
}