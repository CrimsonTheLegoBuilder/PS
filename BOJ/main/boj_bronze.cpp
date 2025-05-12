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

int N, M, K, T;
ll A, B, C;
int I[1005];
std::vector<ll> V;
std::string S;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	while (std::cin >> S) {
		int sz = S.size();
		for (int i = 0; i < sz; i++) {
			if (i == 0 || S[i - 1] != S[i]) std::cout << S[i];
		}
		std::cout << " ";
	} 
	return 0;
}