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
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;

ll p, q;
std::string S, s;
int I[10];
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cin >> p >> q;
	if (p <= 50 && q <= 10) std::cout << "White\n";
	else if (q >= 30) std::cout << "Red\n";
	else std::cout << "Yellow\n";
	return 0;
}