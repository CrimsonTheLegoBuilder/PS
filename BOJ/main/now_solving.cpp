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
const ld PI = acos(-1);

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, I;
int F[8];
char S[15];
//ld p, p0, X, Y, D, A, B, F;
//std::string S;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	N = 5;
	T = 1;
	while (N--) {
		std::cin >> S;
		for (int i = 0; i <= 10; i++) {
			if (S[i + 2] == '\0') break;
			if (S[i] == 'F' && S[i + 1] == 'B' && S[i + 2] == 'I') std::cout << T << " ", M++;
		}
		T++;
	}
	if (!M) std::cout << "HE GOT AWAY!\n";
	return 0;
}