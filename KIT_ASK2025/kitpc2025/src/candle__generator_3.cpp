#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
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

int pt;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
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
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << "(" << p.x << ", " << p.y << ")"; return os; }
	void println() const { std::cout << x << " " << y << "\n";  return; }
	void print() const { std::cout << x << " " << y;  return; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Vpos;
Pos qry[605][3];
Vpos P[605];
bool cmp(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ll tq = p / q;
	return !tq ? p.Euc() <= q.Euc() : tq > 0;
}
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
std::vector<Pos> graham_scan(std::vector<Pos>& C) {
	std::vector<Pos> H;
	if (C.size() < 3) {
		std::sort(C.begin(), C.end());
		return C;
	}
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pos& p, const Pos& q) -> bool {
		int ret = ccw(C[0], p, q);
		if (!ret) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i]) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	return H;
}
int main() {
	int N = 50;
	Pos p = Pos(-1000, 0);
	//std::cout << N << "\n";
	ld t = 2 * PI / N;
	Vpos H;
	//std::cout << N / 3 << "\n";
	for (int i = 0; i < N; i++) {
		Pos p_ = p.rot(t * i);
		//p_.println();
		H.push_back(p_);
	}
	H = graham_scan(H);
	P[0] = H; pt++;
	Pos q = Pos(-900, 0);
	Pos r = Pos(-300, 0);
	ld t_ = t / 6;

	//Vpos T1 = { Pos(-63, -200), Pos(-60, -203), Pos(-57, -200) };
	//Vpos T2 = { Pos(137, 127), Pos(140, 122), Pos(143, 127) };
	//Vpos T3 = { Pos(-143, 127), Pos(-140, 122), Pos(-137, 127) };
	
	//Vpos T1 = { Pos(-63, -300), Pos(-60, -303), Pos(-57, -300) };
	//Vpos T2 = { Pos(137 + 85, 127 + 50), Pos(140 + 85, 122 + 50), Pos(143 + 85, 127 + 50) };
	//Vpos T3 = { Pos(-143 - 85, 127 + 50), Pos(-140 - 85, 122 + 50), Pos(-137 - 85, 127 + 50) };

	//Vpos T1 = { Pos(-63, -400), Pos(-60, -403), Pos(-57, -400) };
	//Vpos T2 = { Pos(137 + 170, 127 + 100), Pos(140 + 170, 122 + 100), Pos(143 + 170, 127 + 100) };
	//Vpos T3 = { Pos(-143 - 170, 127 + 100), Pos(-140 - 170, 122 + 100), Pos(-137 - 170, 127 + 100) };

	//Vpos T1 = { Pos(-63, -500), Pos(-60, -503), Pos(-57, -500) };
	//Vpos T2 = { Pos(137 + 255, 127 + 150), Pos(140 + 255, 122 + 150), Pos(143 + 255, 127 + 150) };
	//Vpos T3 = { Pos(-143 - 255, 127 + 150), Pos(-140 - 255, 122 + 150), Pos(-137 - 255, 127 + 150) };


	//for (int i = 0; i < 11; i++) {
	//	P[pt] = T1; pt++;
	//	P[pt] = T2; pt++;
	//	P[pt] = T3; pt++;
	//	for (Pos& p : T1) p += Pos(12, 0);
	//	for (Pos& p : T2) p += Pos(6, -10);
	//	for (Pos& p : T3) p += Pos(-6, -10);
	//}
	//std::cout << "H = [";
	//for (int i = 0; i < pt; i++) {
	//	std::cout << "  [";
	//	int sz = P[i].size();
	//	for (int j = 0; j < sz; j++) {
	//		std::cout << P[i][j];
	//		if (j < sz - 1) std::cout << ", ";
	//	}
	//	std::cout << "],\n";
	//}
	//std::cout << "]";
	//std::cout << "FUCK::\n";

	//int x_ = -60;
	//Pos p0 = Pos(x_, -100);
	//Pos p1 = Pos(x_ + 1, -101);
	//Pos p2 = Pos(x_ + 2, -100);
	//std::cout << "16\n";
	//for (int i = 1; i <= 16; i++) {
	//	P[i] = { p0, p1, p2 };
	//	p0.x += 6;
	//	p1.x += 6;
	//	p2.x += 6;
	//	//std::cout << "3\n";
	//	//for (Pos& p : P[i]) (p).println();
	//}

	N = 50;
	p = Pos(-500, 0);
	//std::cout << N << "\n";
	t = 2 * PI / N;
	H.clear();
	//std::cout << N / 3 << "\n";
	for (int i = 0; i < N; i++) {
		Pos p_ = p.rot(t * i);
		//p_.println();
		H.push_back(p_);
	}
	P[1] = graham_scan(H);

	freopen("../tests/candle_renew/in/95.in", "w", stdout);

	std::cout << P[0].size() << "\n";
	for (Pos& p : P[0]) p.println();
	std::cout << "1\n";
	std::cout << P[1].size() << "\n";
	for (Pos& p : P[1]) p.println();


	Pos p0, p1, p2;
	p0 = Pos(900, 0);
	p1 = Pos(-900, 0);
	p2 = Pos(0, 900);
	(p0).print();
	std::cout << " ";
	(p1).print();
	std::cout << " ";
	(p2).print();
	std::cout << "\n";
	

	//std::cout << "H = [";
	//for (int i = 0; i < 100; i++) {
	//	Pos p_ = p.rot(t * i);
	//	std::cout << "(" << p_ << "), ";
	//}
	//std::cout << "]";
	return 0;
}