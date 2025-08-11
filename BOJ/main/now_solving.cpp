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

int Q, N, M, Y, K, P, C, L, T, n, m, k, wc, hc, ws, hs;
ll a, d, c, g, b, s, t, x, y, w, l, h, p, q;
int q1, q2, q3, q4, q0;
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	ll vk, jk, vl, jl, vh, dh, jh;
	std::cin >> vk >> jk >> vl >> jl >> vh >> dh >> jh;
	std::cout << (vk * jk + vl * jl) * vh * dh * jh << "\n";
	return 0;
}