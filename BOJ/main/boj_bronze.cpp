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
	for (int i = 0; i < 10; i++) {
		std::cin >> I[i];
	}
	for (int i = 0; i < 10; i++) {
		T = 0;
		for (int j = 0; j < 10; j++) {
			if (i == j) continue;
			T += I[j];
		}
		if (T == I[i]) {
			std::cout << I[i] << "\n";
			break;
		}
	}
	return 0;
}