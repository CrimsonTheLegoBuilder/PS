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

int N, M;
int main() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	for (int _ = 0; _ < 4; _++) std::cin >> N, M += N;
	if (M <= 1500) std::cout << "Yes\n";
	else std::cout << "No\n";
	return 0;
}