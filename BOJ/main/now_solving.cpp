#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <cassert>
#include <unordered_map>
typedef unsigned long long ull;
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
const int LEN = 105;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !sign(x); }
//ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

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
	int x, y, t;
	//ll x, y;
	Pos(int x_ = 0, int y_ = 0, int t_ = 0) : x(x_), y(y_), t(t_) {}
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
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3) / (d2 - d1).mag(); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
int collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
bool between(const Pos& d0, const Pos& d1, const Pos& q) { return sign(dot(d0, d1, q)) < 0 && sign(dot(d1, d0, q)) < 0; }
bool intersect(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2, const int& f = STRONG) {
	bool f1 = ccw(p1, p2, q1) * ccw(p2, p1, q2) > 0;
	bool f2 = ccw(q1, q2, p1) * ccw(q2, q1, p2) > 0;
	if (f == WEAK) return f1 && f2;
	bool f3 = on_seg_strong(p1, p2, q1) ||
		on_seg_strong(p1, p2, q2) ||
		on_seg_strong(q1, q2, p1) ||
		on_seg_strong(q1, q2, p2);
	return (f1 && f2) || f3;
}
bool get_intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2, Pos& t, const int& F = STRONG) {
	if (!ccw(p1, p2, q1, q2)) return 0;
	if (p1.x == p2.x) t.x = p1.x;
	else t.x = p2.x;
	if (p1.y == p2.y) t.y = p1.y;
	else t.y = p2.y;
	if (F == STRONG) return on_seg_strong(p1, p2, t) && on_seg_strong(q1, q2, t);
	return on_seg_strong(q1, q2, t);
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(0, 0), Pos e_ = Pos(0, 0)) : s(s_), e(e_) {}
	bool operator == (const Seg& S) const { return s == S.s && e == S.e; }
	bool operator != (const Seg& S) const { return !(*this == S); }
	//bool operator < (const Seg& S) const { return s == S.s ? e < S.e : s < S.s; }
	bool operator < (const Seg& S) const {
		if (dot(s, e, S.s) >= 0) return 1;
		if (dot(S.s, S.e, s) >= 0) return 0;
		if (s == S.s) return (s - e).Euc() < (s - S.e).Euc();
		return (s - e).Euc() >= (S.s - e).Euc();
	}
	friend std::ostream& operator << (std::ostream& os, const Seg& S) { os << S.s << " " << S.e; return os; }
};
typedef std::vector<Seg> Segs;
Polygon box(const Seg& s) {
	int x1 = std::min(s.s.x, s.e.x);
	int x2 = std::max(s.s.x, s.e.x);
	int y1 = std::min(s.s.y, s.e.y);
	int y2 = std::max(s.s.y, s.e.y);
	return { Pos(x1, y1), Pos(x2, y1), Pos(x2, y2), Pos(x1, y2) };
}
bool floor_check(const Pos& p, const Pos& u, const Pos& d, const Seg& b, Pos& f1, Pos& f2) {
	int up = ccw(p, d, u);
	int c1 = ccw(p, d, b.s);
	int c2 = ccw(p, d, b.e);
	if (c1 && c2) return 0;
	if (up == c1 || up == c2) return 0;
	assert(!c1 || !c2);
	Polygon B = box(b);
	for (int i = 0, j; i < 4; i++) {
		j = (i + 1) % 4;
		if (!ccw(p, d, B[i]) && !ccw(p, d, B[j])) {
			f1 = B[i], f2 = B[j];
			if (dot(p, d, f1, f2) > 0) std::swap(f1, f2);
		}
	}
	return 1;
}
Seg get_floor(const Segs& B, const Pos& p, const Pos& u, const Pos& d) {
	Pos f1, f2;
	Segs V;
	//int sz = B.size();
	//for (int i = 0; i < sz; i++) {
	for (const Seg& b : B) {
		if (floor_check(p, u, d, b, f1, f2)) V.emplace_back(f1, f2);
	}
	std::sort(V.begin(), V.end());
	Pos e;
	bool f = 0;
	for (const Seg& v : V) {
		if (dot(v.s, v.e, p) >= 0) continue;
		if (!f) { e = v.e; f = 1; continue; }
		if (dot(p, e, v.s) > 0) { e = v.e; break; }
		if ((p - e).Euc() < (p - v.e).Euc()) e = v.e;
	}
	return Seg(p, e);
}
bool block_check(const Pos& p, const Pos& u, const Pos& d, const Pos& t, const Seg& b, Pos& e) {
	int up = ccw(p, d, u);
	int c1 = ccw(p, d, b.s);
	int c2 = ccw(p, d, b.e);
	if (c1 * c2 > 0) return 0;
	if (up * c1 <= 0 && up * c2 <= 0) return 0;
	if (dot(p, d, b.s) < 0 || dot(p, d, b.e) < 0) return 0;
	Polygon B = box(b);
	Pos x;
	bool f = 0;
	for (int i = 0, j; i < 4; i++) {
		j = (i + 1) % 4;
		if (get_intersection(p, t, B[i], B[j], x)) {
			if (!f || dot(p, t, x) < dot(p, t, e)) {
				f = 1; e = x;
			}
		}
	}
	return f;
}
bool move(const Segs& B, const Seg& fr, Pos& p, Pos& u, Pos& d, const Pos& s) {
	bool blk = 0, fin = 0;
	Pos t = fr.e, e;
	for (const Seg& b : B) if (block_check(p, u, d, t, b, e)) { blk = 1; t = e; }
	int l = std::max(std::abs(p.x - t.x), std::abs(p.y - t.y));
	if (on_seg_weak(fr.s, fr.e, s)) { t = s; fin = 1; }
	t.t = p.t;
	p = t;
	p.t += l;
	u *= -1;
	std::swap(u, d);
	if (blk) u *= -1, d *= -1;
	return fin;
}
struct Pos3D {
	ll x, y, z, t;
	Pos3D(ll x_ = 0, ll y_ = 0, ll z_ = 0, ll t_ = 0) : x(x_), y(y_), z(z_), t(t_) {}
	bool operator == (const Pos3D& p) const { return x == p.x && y == p.y && z == p.z; }
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	ll operator * (const Pos3D& p) const { return x * p.x + y * p.y + z * p.z; }
	Pos3D operator / (const Pos3D& p) const {
		Pos3D ret;
		ret.x = y * p.z - z * p.y;
		ret.y = z * p.x - x * p.z;
		ret.z = x * p.y - y * p.x;
		return ret;
	}
	Pos3D operator * (const ll& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ll& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ll& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ll& n) { x /= n; y /= n; z /= n; return *this; }
	ll Euc() const { return x * x + y * y + z * z; }
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z; return os; }
}; const Pos3D _1 = Pos3D(1, 1, 1);
typedef std::vector<Pos3D> Polyhedron;
Polyhedron P[LEN];//path
Pos3D P0[LEN];//p0
Pos3D NM[LEN];//norm
ll period[LEN];
ull hash(const Pos3D& q) {
	ull h = 0;
	h |= ull(q.x & 0x1FFFFF);
	h |= ull(q.y & 0x1FFFFF) << 21;
	h |= ull(q.z & 0x1FFFFF) << 42;
	return h;
}
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return !zero(cross(d1, d2, d3).Euc()) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).Euc()) && dot(d1, d3, d2) > 0; }
ll dist(const Pos3D& p, const Pos3D& q) { return std::max({ std::abs(p.x - q.x), std::abs(p.y - q.y), std::abs(p.z - q.z) }); }
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
			std::string s; std::cin >> s;
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
	void get_start_pos(Pos& s, Pos& u, Pos& d) const {
		Pos3D tq = n / v;
		if (tq.x) { s = Pos(p.y, p.z); u = Pos(v.y, v.z); d = Pos(v.y, v.z); }
		if (tq.y) { s = Pos(p.z, p.x); u = Pos(v.z, v.x); d = Pos(v.z, v.x); }
		if (tq.z) { s = Pos(p.x, p.y); u = Pos(v.x, v.y); d = Pos(v.x, v.y); }
	}
	bool slice(const Cube& c, Seg& b) const {
		Pos3D tq = n / v;
		if (tq.x && c.a.x <= p.x && p.x <= c.b.x) { b.s = Pos(c.a.y, c.a.z), b.e = Pos(c.b.y, c.b.z); return 1; }
		if (tq.y && c.a.y <= p.y && p.y <= c.b.y) { b.s = Pos(c.a.z, c.a.x), b.e = Pos(c.b.z, c.b.x); return 1; }
		if (tq.z && c.a.z <= p.z && p.z <= c.b.z) { b.s = Pos(c.a.x, c.a.y), b.e = Pos(c.b.x, c.b.y); return 1; }
		return 0;
	}
} R[LEN];
void pos_push_back(const Pos3D& r, const Pos3D& n, const int& k, const Pos& p) {
	if (n.x) P[k].push_back(Pos3D(r.x, p.x, p.y, p.t));
	if (n.y) P[k].push_back(Pos3D(p.y, r.y, p.x, p.t));
	if (n.z) P[k].push_back(Pos3D(p.x, p.y, r.z, p.t));
	return;
}
void get_path(const int& k) {
	Pos s, u, d;
	const Robot& r0 = R[k];
	r0.get_start_pos(s, u, d);
	Pos3D n = R[k].n / R[k].v;
	Segs B; Seg b;
	for (int i = 0; i < N; i++) if (r0.slice(C[i], b)) B.push_back(b);
	assert(B.size());
	Pos cur = s;
	pos_push_back(r0.p, n, k, cur);
	while (1) {
		Seg fr = get_floor(B, cur, u, d);
		Pos prev = cur;
		bool fin = move(B, fr, cur, u, d, s);
		if (fin) break;
		pos_push_back(r0.p, n, k, cur);
	}
	P0[k] = r0.p;
	NM[k] = r0.n / r0.v;
	assert(!(cur.t % 2));
	period[k] = cur.t >> 1;
	return;
}
ll get_intersection_line(const int& i, const int& j, Pos3D& p0, Pos3D& v) {
	v = NM[i] / NM[j];
	if (v.Euc() == 0) {
		Pos3D q = P0[i] + NM[i];
		ll fc = dot(q, P0[i], P0[j]);
		if (fc) return 0;
		int sz = P[i].size();
		for (int k = 0, l; k < sz; k++) {
			l = (k + 1) % sz;
			if (on_seg_strong(P[i][k], P[i][l], R[j].p)) {
				if (NM[i] == NM[j]) return -1;
				ll dst = dist(P[i][k], R[j].p);
				ll t = P[i][k].t + dst;
				return -t;
			}
		}
		return 0;
	}
	const Pos3D& a = P0[i], & b = P0[j];
	const Pos3D va = NM[i] / v, vb = NM[j] / v;
	Pos a1, a2, b1, b2, x;
	if (v.x) {
		a1 = Pos(a.y, a.z); a2 = a1 + Pos(va.y, va.z);
		b1 = Pos(b.y, b.z); b2 = b1 + Pos(vb.y, vb.z);
	}
	else if (v.y) {
		a1 = Pos(a.z, a.x); a2 = a1 + Pos(va.z, va.x);
		b1 = Pos(b.z, b.x); b2 = b1 + Pos(vb.z, vb.x);
	}
	else if (v.z) {
		a1 = Pos(a.x, a.y); a2 = a1 + Pos(va.x, va.y);
		b1 = Pos(b.x, b.y); b2 = b1 + Pos(vb.x, vb.y);
	}
	get_intersection(a1, a2, b1, b2, x);
	if (v.x) p0 = Pos3D(0, x.x, x.y);
	else if (v.y) p0 = Pos3D(x.y, 0, x.x);
	else if (v.z) p0 = Pos3D(x.x, x.y, 0);
	return 1;
}
int intersection(const Pos3D& p1, const Pos3D& p2, const Pos3D& q1, const Pos3D& q2, const Pos3D& n, Pos3D& q) {
	Pos a1, a2, b1, b2, x;
	if (n.x) { a1 = Pos(p1.y, p1.z); a2 = Pos(p2.y, p2.z); b1 = Pos(q1.y, q1.z); b2 = Pos(q2.y, q2.z); }
	if (n.y) { a1 = Pos(p1.z, p1.x); a2 = Pos(p2.z, p2.x); b1 = Pos(q1.z, q1.x); b2 = Pos(q2.z, q2.x); }
	if (n.z) { a1 = Pos(p1.x, p1.y); a2 = Pos(p2.x, p2.y); b1 = Pos(q1.x, q1.y); b2 = Pos(q2.x, q2.y); }
	if (!get_intersection(a1, a2, b1, b2, x, WEAK)) return -1;
	int t = std::max(std::abs(q1.x - x.x), std::abs(q1.y - x.y));
	if (n.x) q = Pos3D(q1.x, x.x, x.y);
	if (n.y) q = Pos3D(x.y, q1.y, x.x);
	if (n.z) q = Pos3D(x.x, x.y, q1.z);
	t += q1.t;
	assert(!(t % 2));
	return t >> 1;
}
ll collision_time(const int& k, const int& l) {
	Pos3D p0, v;
	ll f1 = get_intersection_line(k, l, p0, v);
	if (!f1) return -1;
	if (f1 < 0) return (-f1 + 1) >> 1;
	std::unordered_map<ull, int> MK, ML;
	std::vector<ull> H;
	Pos3D q;
	int szk = P[k].size();
	for (int i = 0, j; i < szk; i++) {
		j = (i + 1) % szk;
		int t = intersection(p0, v, P[k][i], P[k][j], NM[k], q);
		if (t < 0) continue;
		ull h = hash(q);
		MK[h] = t;
		H.push_back(h);
	}
	int szl = P[l].size();
	for (int i = 0, j; i < szl; i++) {
		j = (i + 1) % szl;
		int t = intersection(p0, v, P[l][i], P[l][j], NM[l], q);
		if (t < 0) continue;
		ull h = hash(q);
		ML[h] = t;
		H.push_back(h);
	}
	std::sort(H.begin(), H.end());
	H.erase(unique(H.begin(), H.end()), H.end());
	ll ans = INF;
	for (const ull& h : H) {
		if (MK.count(h) && ML.count(h)) {
			ll tk = MK[h], pk = period[k];
			ll tl = ML[h], pl = period[l];
			ExtGcd ret = ext_gcd(pk, pl);
			ll diff = tl - tk;
			if (diff % ret.g == 0) {
				ll mod = pl / ret.g;
				ll n = (ret.x % mod) * ((diff / ret.g) % mod) % mod;
				if (n < 0) n += mod;
				ll t_col = tk + n * pk;
				ans = std::min(ans, t_col);
			}
		}
	}
	return ans >= INF ? -1 : ans;
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