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
const ld PI = acos(-1);

int Q, N, M, K, a, b, c, T, W, H, L, s, o;
int q1, q2, q3, q4, q0, A, B, C, D, E, P, S, I;
//int F[2001];
//ld p, p0, X, Y, D, A, B, F;
//std::string S;
void cal_mile(ld d, ld f, ld t, ld& m, ld& mph) {
	d /= 12.;
	d /= 5280.;
	d *= f;
	d *= PI;
	m = d;
	mph = d / t * 3600;
	return;
}
int main() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(2);
	N = 1;
	//ld d = 1;
	//d /= 12.;
	//d /= 5280.;
	//std::cout << d << "\n";
	while (1) {
		ld d, t;
		int f;
		std::cin >> d >> f >> t;
		if (!f) break;
		std::cout << "Trip #" << N << ": ";
		ld m, mph;
		cal_mile(d, f, t, m, mph);
		std::cout << m << " " << mph << "\n";
		N++;
	}
	return 0;
}