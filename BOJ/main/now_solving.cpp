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

int Q, N, M, K;
int q1, q2, q3, q4, q0, n[3];
char B[1001];
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> n[0] >> n[1] >> n[2];
	std::sort(n, n + 3);
	if (n[0] + n[1] == n[2]) std::cout << "1\n";
	else if (n[0] * n[1] == n[2]) std::cout << "2\n";
	else std::cout << "3\n";
	return 0;
}