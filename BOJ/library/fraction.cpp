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
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
bool inner_check(const Polygon& H, const Pos& q) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) if (cross(H[i], H[(i + 1) % sz], q) <= 0) return 0;
	return 1;
}
bool inner_check(Pos p0, Pos p1, Pos p2, const Pos& q) {
	if (cross(p0, p1, p2) < 0) std::swap(p1, p2);
	return inner_check({ p0, p1, p2 }, q);
}
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
};
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
	bool operator <= (const Frac& o) const { return num * o.den <= o.num * den; }
	bool operator == (const Frac& o) const { return num == o.num && den == o.den; }
} _1 = Frac(1), _0 = Frac(0);
struct Range {
	Frac s, e;
	Range(Frac s_ = Frac(-1), Frac e_ = Frac(-1)) : s(s_), e(e_) {}
	bool operator < (const Range& r) const { return s == r.s ? e < r.e : s < r.s; }
};
Frac intersection(const Seg& s1, const Seg& s2) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ll det = (q2 - q1) / (p2 - p1);
	if (!det) return Frac(-1);
	ll a1 = (q2 - q1) / (q1 - p1);
	Frac f1 = Frac(a1, det); f1.fit(); return f1;
	//ll a2 = (p2 - p1) / (p1 - q1);
	//Frac f2 = Frac(-a2, det);
	//if (_0 < f1 && f1 < _1 && _0 <= f2 && f2 <= _1) return f1;
	//return Frac(-1);
}
Range range(const Seg& s1, const Seg& s2, const Pos& o) {
	Pos p0 = o, p1 = s1.s, p2 = s1.e;
	Pos q1 = s2.s, q2 = s2.e;
	if (cross(p0, p1, p2) < 0) std::swap(p1, p2);
	if (cross(p0, q1, q2) < 0) std::swap(q1, q2);
	if (intersect(p0, p1, s2.s, s2.e) && intersect(p0, p2, s2.s, s2.e))
		return Range(_0, _1);
	bool i1 = inner_check(p0, p1, p2, q1);
	bool i2 = inner_check(p0, p1, p2, q2);
	if (!i1 && !i2) return Range();
	Frac f1 = _0, f2 = _1;
	if (i1) f1 = intersection(s1, Seg(p0, q1));
	if (i2) f2 = intersection(s1, Seg(p0, q2));
	return Range(f1, f2);
}