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

int Q, N, M, P, a, b, s, t;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> M >> a >> b;
	a--; b--;
	t = 0;
	while (a != b) {
		a++; t++;
		a = a % M;
	}
	std::cout << t << "\n";
	return 0;
}