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

int Q, N, M, Y, K, P, C, L, T, n, m, k;
ll a, d, c, g, b, s, t, x, y, w, l, h;
int q1, q2, q3, q4, q0;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N >> M >> K;
	T = std::min(M, N / 2);
	//std::cout << "T:: " << T << "\n";
	n = N - T * 2;
	m = M - T;
	//std::cout << "n:: " << n << "\n";
	//std::cout << "m:: " << m << "\n";
	K -= n;
	K -= m;
	K = std::max(K, 0);
	k = K % 3;
	K /= 3;
	//std::cout << "k:: " << k << "\n";
	//std::cout << "K:: " << K << "\n";
	K += k ? 1 : 0;
	T -= K;
	std::cout << T << "\n";
	return 0;
}