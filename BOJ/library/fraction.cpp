#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline ll gcd(ll a, ll b) { while (b) { a %= b; std::swap(a, b); }return a; }

struct Frac {
	ll num, den;//numerator, denominator (num/den)
	Frac(ll n = 0, ll d = 1) : num(n), den(d) {
		assert(den);
		if (den < 0) num *= -1, den *= -1;
		simplify();
	}
	void fit() {
		if (num < 0) num = 0, den = 1;
		else if (num > den) num = den = 1;
		return;
	}
	void simplify() {
		if (!num) { den = 1; return; }
		ll com = gcd(std::abs(num), std::abs(den));
		num /= com; den /= com;
		return;
	}
	friend std::ostream& operator << (std::ostream& os, const Frac& f) {
		os << f.num; if (f.den != 1) os << "/" << f.den; return os;
	}
	bool operator < (const Frac& o) const { return num * o.den < o.num * den; }
	bool operator == (const Frac& o) const { return num == o.num && den == o.den; }
} _1 = Frac(1), _0 = Frac(0);
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
};
struct Range {
	Frac s, e;
	Range(Frac s_ = Frac(), Frac e_ = Frac()) : s(s_), e(e_) {}
	bool operator < (const Range& r) const { return s == r.s ? e < r.e : s < r.s; }
};
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
};
Frac intersection(const Seg& s1, const Seg& s2) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ll det = (q2 - q1) / (p2 - p1);
	if (!det) return Frac(-1);
	//ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ll a1 = (q2 - q1) / (q1 - p1);
	//ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	ll a2 = (p2 - p1) / (p1 - q1);
	//if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return Frac(-1);
}