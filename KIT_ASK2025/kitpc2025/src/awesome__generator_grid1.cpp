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
	N = 2000;
	freopen("../tests/awesome/13.in", "w", stdout);
	std::cout << W << " " << H << " " << N << "\n";
	for (int i = 0; i < 50; i++) {
		int xl = i * 2e7, xr = (i + 1) * 2e7;
		for (int j = 0; j < 40; j++) {
			int yl = j * 25e6, yh = (j + 1) * 25e6;
			
			int x0 = rnd.next(xl, xr - 1);
			int y0 = rnd.next(yl, yh - 1);
			int x1 = rnd.next(x0 + 1, xr);
			int y1 = rnd.next(y0 + 1, yh);
			std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		}
	}
	return 0;
}