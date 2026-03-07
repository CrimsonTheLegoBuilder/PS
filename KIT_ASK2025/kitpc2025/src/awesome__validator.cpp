#define _CRT_SECURE_NO_WARNINGS
#include "../inc/testlib.h"
//#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
typedef long long ll;

int main(int argc, char* argv[]) {
	registerValidation(argc, argv);

	int W, H, N;
	W = inf.readInt(1, 1'000'000'000, "W");
	inf.readSpace();
	H = inf.readInt(1, 1'000'000'000, "H");
	inf.readSpace();
	N = inf.readInt(1, 2000, "N");
	//inf.readEoln();
	inf.skipBlanks();

	int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
	for (int i = 0; i < N; i++) {
		x0 = inf.readInt(0, W - 1, "x0");
		inf.readSpace();
		y0 = inf.readInt(0, H - 1, "y0");
		inf.readSpace();
		x1 = inf.readInt(x0 + 1, W, "x1");
		inf.readSpace();
		y1 = inf.readInt(y0 + 1, H, "y1");
		//inf.readEoln();
		inf.skipBlanks();
	}

	inf.readEof();

	return 0;
}