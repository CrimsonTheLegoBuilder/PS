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

int N, M, K;
ll a, b, C, T;
std::vector<ll> V;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	V.resize(4);
	for (int i = 0; i < 4; i++) std::cin >> V[i];
	std::sort(V.begin(), V.end());
	std::cout << (V[1] + V[2] + V[3] + 1) << "\n";
	return 0;
}