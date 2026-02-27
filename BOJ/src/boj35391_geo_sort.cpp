#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
//ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

int N;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pos& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	Pos operator ~ () const { return { -y, x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
struct Line {
	ll a, b, c;
	bool operator < (const Line& l) const { return a == l.a ? b == l.b ? c < l.c : b < l.b : a < l.a; }
	bool operator == (const Line& l) const { return a == l.a && b == l.b && c == l.c; }
};
Line bisector(const Pos& s, const Pos& e) {
	ll a = 2ll * (e.x - s.x);
	ll b = 2ll * (e.y - s.y);
	ll c = -((ll)e.x * e.x + (ll)e.y * e.y - (ll)s.x * s.x - (ll)s.y * s.y);
	ll g = gcd(std::abs(a), gcd(std::abs(b), std::abs(c)));
	a /= g; b /= g; c /= g;
	if (a < 0 || (a == 0 && b < 0)) a *= -1, b *= -1, c *= -1;
	return { a, b, c };
}
typedef std::vector<Line> Vline;
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(6);
	int R;
	std::cin >> N >> R;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	Vline L;
	for (int i = 0; i < N; i++) 
		for (int j = i + 1; j < N; j++) 
			L.push_back(bisector(P[i], P[j]));
	int cnt = 0;
	int sz = L.size();
	std::sort(L.begin(), L.end());
	for (int i = 0, j; i < sz; i = j) {
		j = i;
		while (j < sz && L[i] == L[j]) j++;
		int n = 2 * (j - i);
		if (n > cnt) cnt = n;
	}
	if (cnt * 100ll >= N * R) std::cout << "YES\n";
	else std::cout << "NO\n";
	return;
}
int main() { solve(); return 0; }//boj35391