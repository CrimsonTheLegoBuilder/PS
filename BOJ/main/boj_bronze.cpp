#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;

int N;
ll X, Y, Z, B, P;
int I[10];
ld a;
std::string S, R;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	std::cin >> a;
	std::cout << 100 / a << "\n" << 100 / (100 - a) << "\n";
	return 0;
}