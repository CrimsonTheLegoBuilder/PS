#define _CRT_SECURE_NO_WARNINGS
#include "../inc/testlib.h"
//#include "testlib.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
typedef long long ll;
typedef long double ld;
typedef std::vector<bool> Vbool;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::vector<ld> Vld;
const ld PI = acos(-1);
inline int sign(const ld& x) { return x < 0 ? -1 : !!x; }
inline ld norm(ld th) {
	while (th < 0) th += 2 * PI;
	while (sign(th - 2 * PI) >= 0) th -= 2 * PI;
	return th;
}
inline int fit(const int& x, const int& lo, const int& hi) { return std::max(lo, std::min(hi, x)); }
inline ll sq(const int& x) { return 1ll * x * x; }

int pt;
int R, L, N, M;
ll D;
struct Pos {
	int x, y, r, w;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) { r = w = 0; }
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	int Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	Pos rot(const ld& the) const { return { int(x * cos(the) - y * sin(the)), int(x * sin(the) + y * cos(the)) }; }
	int quad() const { return y > 0 || y == 0 && x >= 0; }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
	void println() const { std::cout << x << " " << y << "\n";  return; }
	void print() const { std::cout << x << " " << y;  return; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Vpos;
Vpos P, Q;
int C[4] = { 1000, 5000, 10000, 50000 };
int RND;
int main(int argc, char* argv[]) {
	registerGen(argc, argv, 0);
	std::cout << "generating start\n";

	// 이 for문이 파일 생성을 반복합니다.
	for (int i = 83; i <= 90; i++) {
		RND = i;
		for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
		// --- 파일 이름 생성 시작 ---
		std::stringstream ss;
		ss << "../tests/duck/in_re/" << std::setw(2) << std::setfill('0') << i << ".in";
		std::string filename = ss.str();
		// --- 파일 이름 생성 끝 ---

		//std::cout << "Creating file: " << filename << "\n"; // 진행 상황 확인용 출력

		// 생성된 파일 이름으로 출력 방향을 바꿉니다.
		freopen(filename.c_str(), "w", stdout);

		// 기존 데이터 생성 로직 (P와 Q를 매번 새로 생성)
		Vpos P, Q; // 벡터를 루프 안에서 초기화하여 이전 데이터가 남지 않게 함
		R = rnd.next(2, int(1e9));
		L = rnd.next(1, R - 1);
		N = 100;
		M = 100;
		D = sq(R);
		while (P.size() < N) {
			for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
			int x = rnd.next(-int(1e9), int(1e9));
			int y = rnd.next(-int(1e9), int(1e9));
			Pos p = Pos(x, y);
			ll d = p.Euc();
			if (d < D) p.r = rnd.next(1, int(1e9)), P.push_back(p);
		}
		while (Q.size() < M) {
			for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
			int x = rnd.next(-int(1e9), int(1e9));
			int y = rnd.next(-int(1e9), int(1e9));
			Pos q = Pos(x, y);
			ll d = q.Euc();
			if (d < D) q.w = C[rnd.next(0, 3)], Q.push_back(q);
		}

		// 파일에 데이터 출력
		std::cout << R << " " << L << "\n";
		std::cout << N << " " << M << "\n";
		for (int j = 0; j < N; j++) std::cout << P[j] << " " << P[j].r << "\n";
		for (int j = 0; j < M; j++) std::cout << Q[j] << "\n";// << Q[j].w << "\n";

		// fclose(stdout); // stdout을 다시 콘솔로 돌리고 싶을 때 사용
	}

	std::cout << "All files generated.\n"; // 이 메시지는 마지막 파일에 쓰이거나, fclose 후 콘솔에 출력됩니다.
	return 0;
}
//int main(int argc, char* argv[]) {
//	registerGen(argc, argv, 0);
//	std::cout << "generating start\n";
//	for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
//	R = rnd.next(2, int(1e9));
//	//R = int(1e9);
//	L = rnd.next(1, R - 1);
//	//L = R - 1;
//	//N = rnd.next(50, 100);
//	N = 100;
//	//M = rnd.next(50, 100);
//	M = 100;
//	D = sq(R);
//	while (P.size() < N) {
//		for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
//		int x = rnd.next(-int(1e9), int(1e9));
//		int y = rnd.next(-int(1e9), int(1e9));
//		Pos p = Pos(x, y);
//		ll d = p.Euc();
//		if (d < D) p.r = rnd.next(1, R - 1), P.push_back(p);
//	}
//	while (Q.size() < M) {
//		for (int _ = 0; _ < RND; _++) rnd.next(1, 10);
//		int x = rnd.next(-int(1e9), int(1e9));
//		int y = rnd.next(-int(1e9), int(1e9));
//		Pos q = Pos(x, y);
//		ll d = q.Euc();
//		if (d < D) q.w = C[rnd.next(0, 3)], Q.push_back(q);
//	}
//	std::cout << "generating done\n";
//	std::cout << "\n\n\n";
//	freopen("../tests/duck/in/20.in", "w", stdout);
//	std::cout << R << " " << L << "\n";
//	std::cout << N << " " << M << "\n";
//	for (int i = 0; i < N; i++) std::cout << P[i] << " " << P[i].r << "\n";
//	for (int i = 0; i < M; i++) std::cout << Q[i] << " " << Q[i].w << "\n";
//	return 0;
//}