//#define _CRT_SECURE_NO_WARNINGS
//#include "../inc/testlib.h"
////#include "testlib.h"
//#include <iostream>
//#include <algorithm>
//#include <vector>
//#include <cmath>
//typedef long long ll;
//
//int W, H, N;
//#define RND 7
//int main(int argc, char* argv[]) {
//	registerGen(argc, argv, 1);
//	std::cout << "generating start\n\n\n";
//	for (int i = 0; i < RND; i++ ) rnd.next(1, 10);
//	//W = rnd.next(1, int(1e9));
//	W = 1e9;
//	//H = rnd.next(1, int(1e9));
//	H = 1e9;
//	//N = rnd.next(0, 2000);
//	N = 50;
//	freopen("../tests/awesome/in/14.in", "w", stdout);
//	std::cout << W << " " << H << " " << N << "\n";
//	for (int _ = 0; _ < N; _++) {
//		for (int i = 0; i < RND; i++) rnd.next(0, 10);
//		int x0 = rnd.next(0, W - 1);
//		int y0 = rnd.next(0, H - 1);
//		int x1 = rnd.next(x0 + 1, W);
//		int y1 = rnd.next(y0 + 1, H);
//		std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
//	}
//	return 0;
//}

#define _CRT_SECURE_NO_WARNINGS
#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

typedef long long ll;

int W, H, N;
int RND;
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 1);
	std::cout << "generating start\n\n\n";

	int num_files_to_generate = 80;

	for (int i = 73; i <= num_files_to_generate; i++) {
		RND = i;
		for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
		std::stringstream ss;
		ss << "../tests/awesome/in/" << std::setw(2) << std::setfill('0') << i << ".in";
		std::string filename = ss.str();

		freopen(filename.c_str(), "w", stdout);

		W = rnd.next(1, int(1e9));
		H = rnd.next(1, int(1e9));
		N = rnd.next(1500, 2000);

		std::cout << W << " " << H << " " << N << "\n";

		for (int j = 0; j < N; j++) {
			for (int _ = 0; _ < RND; _++) rnd.next(1, 10);

			int x0 = rnd.next(0, W - 1);
			int y0 = rnd.next(0, H - 1);
			int x1 = rnd.next(x0 + 1, W);
			int y1 = rnd.next(y0 + 1, H);
			std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
		}
	}

	//std::cout << "All files generated.\n";
	return 0;
}