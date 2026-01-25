#define _CRT_SECURE_NO_WARNINGS 
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ld INF = 1e17;
const ld TOL = 1e-10;
const ld PI = acos(-1);
const int LEN = 10005;
inline int sign(const ld& x) { return x <= -TOL ? -1 : x >= TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ll sq(const ll& x) { return x * x; }
inline ld sq(const ld& x) { return x * x; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }

#define INSIDE 0
#define MEET 1
#define OUTSIDE 2
#define STRONG 0
#define WEAK 1
#define LO x
#define HI y

int T, N, M;
bool F[3];
Vld Q;
Vint sts;
struct Pos {
	ld x, y;
	Pos(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	bool operator <= (const Pos& p) const { return *this < p || *this == p; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& n) const { return { x * n, y * n }; }
	Pos operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pos& p) const { return x * p.x + y * p.y; }
	ld operator / (const Pos& p) const { return x * p.y - y * p.x; }
	Pos operator ^ (const Pos& p) const { return { x * p.x, y * p.y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ld xy() const { return x * y; }
	Pos rot(const ld& t) const { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
	bool close(const Pos& p) const { return zero((*this - p).Euc()); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
bool inner_check(const Polygon& H, const Pos& p) {
	int sz = H.size();
	for (int i = 0; i < sz; i++) if (ccw(H[i], H[(i + 1) % sz], p) < 0) return 0;
	return 1;
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(), Pos e_ = Pos()) : s(s_), e(e_) {}
	Pos p(const ld& rt) const { return s + (e - s) * rt; }
	ld green(const ld& lo = 0, const ld& hi = 1) const {
		ld d = hi - lo;
		ld ratio = (lo + hi) * .5;
		Pos m = p(ratio);
		return m.y * d * (s.x - e.x);
	}
};
ld intersection(const Seg& s1, const Seg& s2, const bool& f = STRONG) {
	const Pos& p1 = s1.s, p2 = s1.e, q1 = s2.s, q2 = s2.e;
	ld det = (q2 - q1) / (p2 - p1);
	if (zero(det)) return -1;
	ld a1 = ((q2 - q1) / (q1 - p1)) / det;
	ld a2 = ((p2 - p1) / (p1 - q1)) / -det;
	if (f == WEAK) return fit(a1, 0, 1);
	if (0 < a1 && a1 < 1 && -TOL < a2 && a2 < 1 + TOL) return a1;
	return -1;
}
struct Circle {
	Pos c;
	ll r;
	Circle(Pos c_ = Pos(), ll r_ = 0) : c(c_), r(r_) {}
	bool operator > (const Pos& p) const { return sign(r - (c - p).mag()) > 0; }
	Pos p(const ld& t) const { return c + Pos(r, 0).rot(t); }
	ld rad(const Pos& p) const { return (p - c).rad(); }
	ld area(const ld& lo, const ld& hi) const { return (hi - lo) * r * r * .5; }
	ld green(const ld& lo, const ld& hi) const {
		Pos s = Pos(cos(lo), sin(lo)), e = Pos(cos(hi), sin(hi));
		ld fan = area(lo, hi);
		Pos m = c + (s + e) * r * (ld).5;
		ld tz = (cos(lo) - cos(hi)) * m.y * r;
		return fan + tz - (s / e) * r * r * (ld).5;
	}
	friend std::istream& operator >> (std::istream& is, Circle& p) { is >> p.c.x >> p.c.y >> p.r; return is; }
	friend std::ostream& operator << (std::ostream& os, const Circle& p) { os << p.c.x << " " << p.c.y << " " << p.r; return os; }
} C[3];
Vld intersections(const Circle& a, const Circle& b) {
	Pos ca = a.c, cb = b.c;
	Pos vec = cb - ca;
	ld ra = a.r, rb = b.r;
	ld distance = vec.mag();
	ld rd = vec.rad();
	if (vec.Euc() > sq(ra + rb) + TOL) return {};
	if (vec.Euc() < sq(ra - rb) - TOL) return {};
	ld X = (ra * ra - rb * rb + vec.Euc()) / (2 * distance * ra);
	if (X < -1) X = -1;
	if (X > 1) X = 1;
	ld h = acos(X);
	Vld ret = {};
	ret.push_back(norm(rd - h));
	if (zero(h)) return ret;
	ret.push_back(norm(rd + h));
	return ret;
}
struct Sphere {
	ll x, y, z, r;
	Sphere(ll x_ = 0, ll y_ = 0, ll z_ = 0, ll r_ = 0) : x(x_), y(y_), z(z_), r(r_) {}
	//bool operator < (const Sphere& q) const { return r > q.r; }
	bool operator < (const Sphere& q) const { return x == q.x ? y == q.y ? z < q.z : y < q.y : x < q.x; }
	Sphere operator - (const Sphere& q) const { return { x - q.x, y - q.y, z - q.z, 0 }; }
	ll operator * (const Sphere& q) const { return x * q.x + y * q.y + z * q.z; }
	ld vol() const { return (4. / 3) * PI * r * r * r; }
	ld vol(const ld& h) const { return PI * h * h * (3 * r - h) / 3; }
	ld surf() const { return PI * 4 * r * r; }
	ld surf(const ld& h) const { return PI * 2 * r * h; }
} S[3], SP[LEN * 3];
bool collinear() {
	ll x1 = S[1].x - S[0].x;
	ll y1 = S[1].y - S[0].y;
	ll z1 = S[1].z - S[0].z;
	ll x2 = S[2].x - S[1].x;
	ll y2 = S[2].y - S[1].y;
	ll z2 = S[2].z - S[1].z;
	return (y1 * z2 - z1 * y2) == 0 && (z1 * x2 - x1 * z2) == 0 && (x1 * y2 - y1 * x2) == 0;
}
ll Euc(const Sphere& p, const Sphere& q) { return sq(p.x - q.x) + sq(p.y - q.y) + sq(p.z - q.z); }
ld mag(const Sphere& p, const Sphere& q) { return sqrtl(Euc(p, q)); }
ld rad(const Sphere& a, const Sphere& b, const Sphere& c) {
	ld dab = mag(a, b);
	ld dac = mag(a, c);
	ll det = (b - a) * (c - a);
	ld proj = det / dab;
	ld ret = fit(proj / dac, -1, 1);
	return acos(ret);
}
int meet(const Sphere& p, const Sphere& q) {
	ll dist = Euc(p, q);
	ll rout = sq(p.r + q.r);
	if (dist >= rout) return OUTSIDE;
	ll rin = sq(p.r - q.r);
	if (dist <= rin) return INSIDE;
	return MEET;
}
ld cos_2nd(const ld& a, const ld& b, const ld& c) {
	ld num = sq(a) + sq(b) - sq(c);
	ld den = 2 * a * b;
	ld t = num / den;
	return std::abs(acosl(std::min(std::max(t, -(ld)1.0), (ld)1.0)));
}
void spherical_triangle_angles(const ld& a, const ld& b, const ld& c, ld& A_, ld& B_, ld& C_) {
	A_ = acos(fit((cos(a) - cos(b) * cos(c)) / (sin(b) * sin(c)), -1, 1));
	B_ = acos(fit((cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c)), -1, 1));
	C_ = acos(fit((cos(c) - cos(a) * cos(b)) / (sin(a) * sin(b)), -1, 1));
	return;
}
ld area(const ld& a, const ld& b, const ld& c, const ll& r, const ld& t) {
	ld A_, B_, C_;
	if (a >= PI) {
		spherical_triangle_angles(a * .5, t, c, A_, B_, C_);
		return r * r * (A_ + B_ + C_ - PI) * 2;
	}
	spherical_triangle_angles(a, b, c, A_, B_, C_);
	return r * r * (A_ + B_ + C_ - PI);
}
bool inner_check(const Circle& c, const Pos& p1, const Pos& p2, const Pos q1, const Pos& q2) {

}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	for (int i = 0; i < 3; i++) std::cin >> S[i].x >> S[i].y >> S[i].z >> S[i].r;
	std::sort(S, S + 3);
	if (collinear()) {
		//평행한 경우 따로 처리
		//한 원이 먹혀있거나 두 고깔 모양
		return;
	}
	//접평면 검사 방법이 필요함
	//그냥 고깔 2개 이어붙인 듯한 모양으로도 만들어질 수 있는 거라서
	//접선 2개를 구한 후 사이에 낀 원이 두 선분에 모두 걸치면 해당 원을 가운데에 둔 두 고깔 형태의 볼록 껍질임\
	//한 원이 먹혀있는지 검사
	for (int i = 0, j, k; i < 3; i++) {
		j = (i + 1) % 3;
		k = (j + 1) % 3;
		ld dab = mag(S[i], S[j]);
		ld dac = mag(S[i], S[k]);
		Pos ca = Pos(0, 0);
		Pos cb = Pos(dab, 0);
		ld t = rad(S[i], S[j], S[k]);
		Pos cc = Pos(dac, 0).rot(t);
		C[0] = Circle(ca, S[i].r);
		C[1] = Circle(cb, S[j].r);
		C[2] = Circle(cc, S[k].r);
		ld dif = S[j].r - S[i].r;
		bool f = 0;
		if (sign(dif) < 0) f = 1, dif *= -1;
		ld tt = acos(f / dab);
		if (f) tt = PI - tt;
		Pos p1 = Pos(S[i].r, 0).rot(tt);
		Pos p2 = cb + p1.unit() * S[j].r;
	}
	//반례를 모두 처리했으므로 계산 시작
	//일단 가운데 잘린 삼각기둥은 변의 길이를 모두 구할 수 있으므로 공식이 따로 존재함
	//한 원을 잡고 다른 원 2개를 그려줌
	//접선 2개씩을 그려서 접점 2개를 이어 직선을 그린다
	//두 직선의 교점 == 법선 벡터가 구에 사영된 점
	//잘린 구는 표면의 구면 직사각형의 넓이를 구해서 구해준다
	//구면직사각형을 바로 구하는 건 아니고,
	//원뿔대에 의해 가려진 부분의 부피는 상대적으로 구하기 쉬우므로
	//360도 전체의 가려진 부분을 구하고 나서
	//안쪽은 원뿔을 파주고
	//각도 비례식 세워서 잘려야 하는 부분을 날려줌
	//근데 결국 3개의 원뿔대의 회전 각도는 잘 구하긴 해야함. 위 방법으로 구하면 되겠지 뭐
	//잘린 회전체의 부피는 회전 각도와 회전 단면의 무게 중심을 구할 수 있으므로 바로 구해짐
	return;
}
int main() { solve(); return 0; }//boj23590
