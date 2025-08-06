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
ll a, d, c, g, b, s, t, x, y, w, l, h;
int q1, q2, q3, q4, q0;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> w >> l >> h;
	if (w < l) std::swap(w, l);
	bool f = 1;
	if (l < h + h) f = 0;
	if (w > l + l) f = 0;
	std::cout << (f ? "good\n" : "bad\n");
	return 0;
}