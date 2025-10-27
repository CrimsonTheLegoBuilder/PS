#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
typedef long long ll;
typedef double ld;
//typedef long double ld;
const ld TOL = 1e-7;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }

//#include <unistd.h>
//int iptr = 0, left = 0;
//char ibuf[1 << 20];
//int read() {
//	int x = 0, s = 1;
//	bool started = false;
//	for (;;) {
//		if (iptr >= left) {
//			left = read(0, ibuf, sizeof ibuf);
//			iptr = 0;
//		}
//		if (!left) break;
//		char c = ibuf[iptr++];
//		if (c == '-') {
//			s = -1;
//		}
//		else if (c >= '0' && c <= '9') {
//			x = x * 10 + (c - '0');
//			started = true;
//		}
//		else if (started) {
//			break;
//		}
//	}
//	return s > 0 ? x : -x;
//}

int N;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return sqrtl(Euc()); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
};
const Pos O = Pos(0, 0);
const Pos pp = Pos(1, 1);
const Pos np = Pos(-1, 1);
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ld intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) {
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	return ((q2 - q1) / (q1 - p1)) / det;
}
Pos seg[8][2];
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(19);
	//std::cin >> N; Polygon P(N); for (Pos& p : P) std::cin >> p;
	////std::cin >> N; Polygon P(N); for (Pos& p : P) p.x = read(), p.y = read();
	//Pos r = P[0], ur = P[0], u = P[0], ul = P[0];
	//Pos l = P[0], dl = P[0], d = P[0], dr = P[0];
	//for (int i = 0; i < N; i++) {
	//	if (u.y < P[i].y) u = P[i];
	//	if (d.y > P[i].y) d = P[i];
	//	if (r.x < P[i].x) r = P[i];
	//	if (l.x > P[i].x) l = P[i];
	//	if (cross(O, pp, ul) < cross(O, pp, P[i])) ul = P[i];
	//	if (cross(O, pp, dr) > cross(O, pp, P[i])) dr = P[i];
	//	if (cross(O, np, dl) < cross(O, np, P[i])) dl = P[i];
	//	if (cross(O, np, ur) > cross(O, np, P[i])) ur = P[i];
	//}
	std::cin >> N; Pos p; std::cin >> p;
	//std::cin >> N; Polygon P(N); for (Pos& p : P) p.x = read(), p.y = read();
	Pos r = p, ur = p, u = p, ul = p;
	Pos l = p, dl = p, d = p, dr = p;
	for (int i = 0; i < N; i++) {
		std::cin >> p;
		if (u.y < p.y) u = p;
		if (d.y > p.y) d = p;
		if (r.x < p.x) r = p;
		if (l.x > p.x) l = p;
		if (cross(O, pp, ul) < cross(O, pp, p)) ul = p;
		if (cross(O, pp, dr) > cross(O, pp, p)) dr = p;
		if (cross(O, np, dl) < cross(O, np, p)) dl = p;
		if (cross(O, np, ur) > cross(O, np, p)) ur = p;
	}
	r = r + Pos(1, 0);
	ur = ur + Pos(1, 0);
	u = u + Pos(0, 1);
	ul = ul + Pos(0, 1);
	l = l + Pos(-1, 0);
	dl = dl + Pos(-1, 0);
	d = d + Pos(0, -1);
	dr = dr + Pos(0, -1);
	seg[0][0] = r; seg[0][1] = r + Pos(0, 1);
	seg[1][0] = ur; seg[1][1] = ur + Pos(-1, 1);
	seg[2][0] = u; seg[2][1] = u + Pos(-1, 0);
	seg[3][0] = ul; seg[3][1] = ul + Pos(-1, -1);
	seg[4][0] = l; seg[4][1] = l + Pos(0, -1);
	seg[5][0] = dl; seg[5][1] = dl + Pos(1, -1);
	seg[6][0] = d; seg[6][1] = d + Pos(1, 0);
	seg[7][0] = dr; seg[7][1] = dr + Pos(1, 1);
	ld ppp = 0;
	for (int i = 0; i < 8; i++) {
		int i0 = (i + 7) % 8, i1 = i, i2 = (i + 1) % 8;
		ld s = intersection(seg[i1][0], seg[i1][1], seg[i0][0], seg[i0][1]);
		ld e = intersection(seg[i1][0], seg[i1][1], seg[i2][0], seg[i2][1]);
		ppp += (seg[i1][1] - seg[i1][0]).mag() * (e - s);
	}
	std::cout << ppp << "\n";
	//printf("%13lf\n", ppp);
	return;
}
int main() { solve(); return 0; }
