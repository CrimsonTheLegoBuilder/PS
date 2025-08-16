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
int q1, q2, q3, q4, q0, n[3], g, p ,t, A, B, C, D, P, I;
int ret[2001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> A >> B >> C >> D >> P;
	a = A * P;
	b = B + (std::max(0, P - C)) * D;
	std::cout << std::min(a, b) << "\n";
	return 0;
}