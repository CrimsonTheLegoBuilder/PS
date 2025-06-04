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

int T, N;
ll X, Y, A, B, C, D, E, Z, P, a, c, e;
int I[10];
std::string S, R;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	std::cin >> A >> C >> E;
	std::cin >> a >> c >> e;
	if (e < E) std::cout << "E\n";
	else {
		if (c < C) {
			if (C / 2 < c) std::cout << "D\n";
			else std::cout << "E\n";
		}
		else {
			if (a < A) {
				if (A / 2 < a) std::cout << "B\n";
				else std::cout << "C\n";
			}
			else std::cout << "A\n";
		}
	}
	return 0;
}