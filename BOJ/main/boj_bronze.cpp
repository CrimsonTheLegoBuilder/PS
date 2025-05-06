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
ll c, s, i, a;
std::vector<ll> V;
char t;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	std::cin >> N;
	while (N--) {
		std::cin >> t;
		if (t == 'C') c++;
		if (t == 'S') s++;
		if (t == 'I') i++;
		if (t == 'A') a++;
	}
	std::cin >> t;
	if (t == 'C') K = c;
	if (t == 'S') K = s;
	if (t == 'I') K = i;
	if (t == 'A') K = a;
	std::cout << K << "\n";
	return 0;
}