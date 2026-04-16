#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
const int LEN = 105;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
//ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

struct ExtGcd {
	ll g, x, y;
	ExtGcd(ll g_ = 0, ll x_ = 0, ll y_ = 0) : g(g_), x(x_), y(y_) {}
};
ExtGcd ext_gcd(ll a, ll b) {
	ll x0 = 1, x1 = 0, y0 = 0, y1 = 1;
	while (b) {
		ll q = a / b; a %= b; std::swap(a, b);
		x0 -= q * x1; std::swap(x0, x1);
		y0 -= q * y1; std::swap(y0, y1);
	}
	return { a, x0, y0 };
}
ll ans;
int N, K;
struct Pos {
	int x, y;
	//ll x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	//Pos(ll x_ = 0, ll y_ = 0) : x(x_), y(y_) {}
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
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(0, 0), Pos e_ = Pos(0, 0)) : s(s_), e(e_) {}
	bool operator == (const Seg& S) const { return s == S.s && e == S.e; }
	bool operator != (const Seg& S) const { return !(*this == S); }
	bool operator < (const Seg& S) const { return s == S.s ? e < S.e : s < S.s; }
	friend std::ostream& operator << (std::ostream& os, const Seg& S) { os << S.s << " " << S.e; return os; }
};
struct Pos3D {
	ll x, y, z, t;
	Pos3D(ll x_ = 0, ll y_ = 0, ll z_ = 0, ll t_ = 0) : x(x_), y(y_), z(z_), t(t_) {}
	bool operator == (const Pos3D& p) const { return x == p.x && y == p.y && z == p.z; }
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	Pos3D operator * (const ll& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ll& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ll& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ll& n) { x /= n; y /= n; z /= n; return *this; }
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z; return os; }
}; const Pos3D _1 = Pos3D(1, 1, 1);
typedef std::vector<Pos3D> Polyhedron;
Polyhedron P[LEN];//path
struct Cube {
	Pos3D a, b;
	Cube(Pos3D a_ = Pos3D(), Pos3D b_ = Pos3D()) : a(a_), b(b_) {}
	Cube& operator *= (const ll& n) { a *= n; b *= n; return *this; }
	void init() { std::cin >> *this; *this *= 2; }
	friend std::istream& operator >> (std::istream& is, Cube& c) { is >> c.a >> c.b; return is; }
	friend std::ostream& operator << (std::ostream& os, const Cube& c) { os << c.a << " " << c.b; return os; }
} C[LEN];
struct Robot {
	Pos3D p, n, v;
	Robot(Pos3D p_ = Pos3D(), Pos3D n_ = Pos3D(), Pos3D v_ = Pos3D()) : p(p_), n(n_), v(v_) {}
	void init() {
		std::cin >> p; p = p * 2 + _1;
		for (Pos3D* q : { &n, &v }) {
			std::string s;
			std::cin >> s;
			if (s == "x+") *q = Pos3D(1);
			else if (s == "x-") *q = Pos3D(-1);
			else if (s == "y+") *q = Pos3D(0, 1);
			else if (s == "y-") *q = Pos3D(0, -1);
			else if (s == "z+") *q = Pos3D(0, 0, 1);
			else if (s == "z-") *q = Pos3D(0, 0, -1);
		}
		p += n;
		return;
	}
} R[LEN];
void get_path(const int& k) {
	std::vector<Seg> B;
	const Robot& r0 = R[k];
	bool x = r0.n.x || r0.v.x;
	bool y = r0.n.y || r0.v.y;
	bool z = r0.n.z || r0.v.z;
	int ref = !x ? r0.p.x : !y ? r0.p.y : r0.p.z;
	
	for (int i = 0; i < N; i++) {

	}
}
ll collision_time(const int& k, const int& l) {
	//
	return -1;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(6);
	std::cin >> N >> K; ans = -1;
	for (int i = 0; i < N; i++) C[i].init();
	for (int k = 0; k < K; k++) R[k].init();
	for (int k = 0; k < K; k++) get_path(k);
	for (int k = 0; k < K; k++) 
		for (int l = k + 1; l < K; l++) 
			ans = std::max(ans, collision_time(k, l));
	if (ans < 0) std::cout << "ok\n";
	else std::cout << ans << "\n";
	return;
}
int main() { solve(); return 0; }//boj23202