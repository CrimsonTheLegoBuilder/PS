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
	//W = 1000000000;
	//H = 1000000000;
	//N = 2000;
	////TX = 9.5e6;
	////TY = 1.2e7;
	//freopen("../tests/awesome/in/60.in", "w", stdout);
	//std::cout << W << " " << H << " " << N << "\n";
	//for (int i = 0; i < 500; i++) {
	//	std::cout << "0 0 500000000 500000000\n";
	//}
	//for (int i = 0; i < 500; i++) {
	//	std::cout << "0 500000000 500000000 1000000000\n";
	//}
	//for (int i = 0; i < 500; i++) {
	//	std::cout << "500000000 0 1000000000 500000000\n";
	//}
	//for (int i = 0; i < 500; i++) {
	//	std::cout << "500000000 500000000 1000000000 1000000000\n";
	//}

	W = 2001000;
	H = 2000;
	N = 2000;
	//TX = 9.5e6;
	//TY = 1.2e7;
	//freopen("../tests/awesome/in/85.in", "w", stdout);

	//std::cout << W << " " << H << " " << N << "\n";
	int x0 = 0, x1 = 2000;
	int y0 = 1999, y1 = 3999;
	for (int i = 0; i < 2000; i++) {
		std::cout << x0 << " " << y0 << " ";
		std::cout << x1 << " " << y1 << "\n";
		x0++;
		y0--;
		x1++;
		y1--;
	}

	//int x0 = 0, x1 = 1;
	//int y0 = 0, y1 = 1;
	//for (int i = 0; i < 2000; i++) {
	//	std::cout << x0 << " " << y0 << " ";
	//	std::cout << x1 << " " << y1 << "\n";
	//	//x0--;
	//	x1++;
	//	//y0--;
	//	y1++;
	//}

	//for (int i = 0; i < N; i++) {
	//	std::cout << "0 0 1000000000 1000000000\n";
	//}
	//for (int i = 0; i < 50; i++) {
	//	int xl = i * 2e7, xr = (i + 1) * 2e7;
	//	for (int j = 0; j < 40; j++) {
	//		int yl = j * 25e6, yh = (j + 1) * 25e6;

	//		int x0 = xl + TX;
	//		int y0 = yl + TY;
	//		int x1 = xr - TX;
	//		int y1 = yh - TY;
	//		std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
	//	}
	//}
	return 0;
}