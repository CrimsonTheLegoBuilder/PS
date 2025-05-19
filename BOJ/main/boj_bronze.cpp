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

ll T, L, N, H, C, A, B, D, a, b;
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
	std::cin >> A >> B;
	D = 1;
	while (A) {
		a = A % 10;
		b = B % 10;
		T += std::max(a, b) * D;
		D *= 10;
		A /= 10;
		B /= 10;
	}
	std::cout << T << "\n";
	return 0;
}