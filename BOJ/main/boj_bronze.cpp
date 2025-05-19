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

int T, L, N, H, C, A, B, D;
std::string S;
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cin >> T;
	//while (T--) {
	//	std::cin >> N;
	//	while (N--) {
	//		std::cin >> A >> B;
	//		std::cout << A + B << " " << A * B << "\n";
	//	}
	//}
	for (int i = 1; i <= T; i++) {
		std::cin >> N;
		std::cout << "Case #" << i << ": ";
		if (N <= 25) std::cout << "World Finals\n";
		else {
			std::cout << "Round ";
			if (N <= 1000) std::cout << "3\n";
			else if (N <= 4500) std::cout << "2\n";
			else std::cout << "1\n";
		}
	}
	return 0;
}