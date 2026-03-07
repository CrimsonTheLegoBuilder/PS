#define _CRT_SECURE_NO_WARNINGS
//#include "../inc/testlib.h"
#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
typedef long long ll;

int W, H, N;
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 1);
	std::cout << "generating start\n\n\n";
	W = 1e9;
	H = 1e9;
	N = 200;
	int X = 20, Y = 10;
	freopen("../tests/awesome/15.in", "w", stdout);
	std::cout << W << " " << H << " " << N << "\n";
	for (int i = 0; i < X; i++) {
		int xl = i * (1e9 / X), xr = (i + 1) * (1e9 / X);
		for (int j = 0; j < Y; j++) {
			int yl = j * (1e9 / Y), yh = (j + 1) * (1e9 / Y);

			int x0 = rnd.next(xl, xr - 1);
			int y0 = rnd.next(yl, yh - 1);
			int x1 = rnd.next(x0 + 1, xr);
			int y1 = rnd.next(y0 + 1, yh);
			std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		}
	}
	return 0;
}