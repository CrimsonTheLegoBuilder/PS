#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <iomanip>
#include <chrono>

long long N, M, K, L, R, X, Y, W;
std::string Dir = "../tests/duck/";
#include <fstream>
void conv() {
	std::cin >> R >> L;
	//std::cout << R << " " << L << "\n";
	std::cin >> N >> M;
	//std::cout << N << " " << M << "\n";
	for (int i = 0; i < N; i++) {
		std::cin >> X >> Y >> R;
		//std::cout << X << " " << Y << " " << R << "\n";
	}
	int t = 0;
	for (int i = 0; i < M; i++) {
		std::cin >> X >> Y;
		if (R * R == (X * X) + (Y * Y)) t++;
	}
	std::cout << t << "\n";
	return;
}
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	//std::cout << std::scientific;
	std::cout.precision(13);
	for (int i = 1; i <= 82; i++) {
		std::string Qin = Dir + "in/" + (i < 10 ? "0" : "") + std::to_string(i) + ".in";
		std::string Qout = Dir + "val/" + (i < 10 ? "0" : "") + std::to_string(i) + ".in";
		std::ifstream fin(Qin);
		std::ofstream fout(Qout);
		if (!fin || !fout) {
			std::cerr << "open error: " << Qin << " or " << Qout << "\n";
			continue;
		}
		std::streambuf* orig_cin = std::cin.rdbuf();
		std::streambuf* orig_cout = std::cout.rdbuf();
		std::cin.rdbuf(fin.rdbuf());
		std::cout.rdbuf(fout.rdbuf());
		conv();
		std::cin.rdbuf(orig_cin);
		std::cout.rdbuf(orig_cout);
	}
	return 0;
}