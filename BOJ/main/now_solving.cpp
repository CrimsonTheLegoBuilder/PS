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

int Q, N, M, K, a, b, c, T, W, H, s;
int q1, q2, q3, q4, q0, n[3], g, p ,t, A, I;
char B[1001];
std::string S;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N;
	while (N--) {
		int r, e, c;
		std::cin >> r >> e >> c;
		a = e - c;
		b = r;
		if (a > b) std::cout << "advertise\n";
		else if (a < b) std::cout << "do not advertise\n";
		else std::cout << "does not matter\n";
	}
	return 0;
}