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
Pos qry[505][3];
Vpos P[500];
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
	int N = 100;
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
	P[0] = H; pt++;
	Pos q = Pos(-900, 0);
	Pos r = Pos(-300, 0);
	ld t_ = t / 6;
	for (int i = 0; i < N; i++) {
		//std::cout << "3\n";
		if (i % 3 == 0) {
			Pos p0 = q.rot(t * i + t_);
			Pos p1 = q.rot(t * i + t_ * 5);
			Pos p2 = r.rot(t * i + t_ * 3);
			//p0.println();
			//p1.println();
			//p2.println();
			//std::cout << "pt:: " << pt << "\n";
			P[pt] = { p0, p1, p2 }; pt++;
		}
	}
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
	std::cout << N << "\n";
	for (Pos& p : P[0]) p.println();
	std::cout << N / 3 << "\n";
	//std::cout << "FUCK\n";
	for (int i = 1; i <= (N / 3); i++) {
		std::cout << P[i].size() << "\n";
		for (Pos& p : P[i]) p.println();
	}
	Pos R = Pos(-204, 200), G = Pos(-200, 200), B = Pos(-196, 200);
	int Q = 100;
	int w = 8;
	int h = -10;
	std::cout << Q << "\n";
	for (int i = 1; i <= Q; i++) {
		R.print();
		std::cout << " ";
		G.print();
		std::cout << " ";
		B.println();
		R.x += w;
		G.x += w;
		B.x += w;
		if (i % 50 == 0) {
			w *= -1;
			R.y += h;
			G.y += h;
			B.y += h;
		}
	}

	//std::cout << "H = [";
	//for (int i = 0; i < 100; i++) {
	//	Pos p_ = p.rot(t * i);
	//	std::cout << "(" << p_ << "), ";
	//}
	//std::cout << "]";
	return 0;
}