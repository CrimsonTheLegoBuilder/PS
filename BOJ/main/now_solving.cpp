#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <deque>
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
const int LEN = 1e6 + 1;
const ld TOL = 1e-7;
bool zero(const ld& x) { return std::abs(x) < TOL; }
int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
int sign(const ll& x) { return x < 0 ? -1 : !!x; }

int N, M, Q, T;
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(9);

	return;
}
int main() { solve(); return 0; }
