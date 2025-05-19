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

ll P, T, L, N, H, C, A, B, D, a, b;
std::string S;
int I[10];
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	//std::cin >> T;
	//while (T--) {
	//	std::cin >> N;
	//	while (N--) {
	//		std::cin >> A >> B;
	//		std::cout << A + B << " " << A * B << "\n";
	//	}
	//}
	std::cin >> N;
	while (N--) {
		std::cin >> P >> T;
		std::cout << P - (T / 7) + (T / 4) << "\n";
	}
	return 0;
}