#define _CRT_SECURE_NO_WARNINGS
//#include "../inc/testlib.h"
#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
typedef long long ll;

int W, H, N, T;
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 1);
	std::cout << "generating start\n\n\n";
	W = 1e9;
	H = 1e9;
	N = 2000;
	T = 1e6;
	freopen("../tests/awesome/in/36.in", "w", stdout);
	std::cout << W << " " << H << " " << N << "\n";
	int x0 = 0;
	int y0 = 0;
	int x1 = 1e6;
	int y1 = 1e6;
	for (int i = 0; i < (N >> 1) ; i++) {
		std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		x0 += T;
		y0 += T;
		x1 += T;
		y1 += T;
	}
	x0 = 0;
	y0 = 1e9 - T;
	x1 = T;
	y1 = 1e9;
	for (int i = 0; i < (N >> 1); i++) {
		std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		x0 += T;
		y0 -= T;
		x1 += T;
		y1 -= T;
	}
	return 0;
}