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

int Q, N, M, K, P, C, L, T, a, d, g, b, s, t;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N >> M >> K;
	if (N + M == K) { std::cout << N << "+" << M << "=" << K << "\n"; return 0; }
	if (N - M == K) { std::cout << N << "-" << M << "=" << K << "\n"; return 0; }
	if (N * M == K) { std::cout << N << "*" << M << "=" << K << "\n"; return 0; }
	if (N / M == K && !(N % M)) { std::cout << N << "/" << M << "=" << K << "\n"; return 0; }
	if (N == M + K) { std::cout << N << "=" << M << "+" << K << "\n"; return 0; }
	if (N == M - K) { std::cout << N << "=" << M << "-" << K << "\n"; return 0; }
	if (N == M * K) { std::cout << N << "=" << M << "*" << K << "\n"; return 0; }
	//if (N == M / K && !(M % K)) { std::cout << N << "=" << M << "/" << K << "\n"; return 0; }
	return 0;
}