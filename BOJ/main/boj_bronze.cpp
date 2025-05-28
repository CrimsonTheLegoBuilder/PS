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

int N;
ll T, S, M, Ds, Ys, Dm, Ym;
int I[10];
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cin >> Ds >> Ys >> Dm >> Ym;
	S = -Ds;
	while (S < 0) S += Ys;
	M = -Dm;
	while (M < 0) M += Ym;
	while (1) {
		while (M < S) M += Ym;
		if (M == S) {
			T = S;
			break;
		}
		S += Ys;
	}
	std::cout << T << "\n";
	return 0;
}