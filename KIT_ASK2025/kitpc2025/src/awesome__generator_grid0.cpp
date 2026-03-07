#define _CRT_SECURE_NO_WARNINGS
//#include "../inc/testlib.h"
#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
typedef long long ll;

int W, H, N, TX, TY;
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 1);
	std::cout << "generating start\n\n\n";
	W = 1e9;
	H = 1e9;
	N = 2000;
	TX = 9.5e6;
	TY = 1.2e7;
	freopen("../tests/awesome/in/30.in", "w", stdout);
	std::cout << W << " " << H << " " << N << "\n";
	for (int i = 0; i < 50; i++) {
		int xl = i * 2e7, xr = (i + 1) * 2e7;
		for (int j = 0; j < 40; j++) {
			int yl = j * 25e6, yh = (j + 1) * 25e6;

			int x0 = xl + TX;
			int y0 = yl + TY;
			int x1 = xr - TX;
			int y1 = yh - TY;
			std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		}
	}
	return 0;
}