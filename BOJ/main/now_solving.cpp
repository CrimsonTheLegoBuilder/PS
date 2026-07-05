#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
typedef unsigned long long ull;
typedef long long ll;
typedef double ld;
const ll INF = 1e17;
const int LEN = 105;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !sign(x); }

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
char dir[5];
struct Pos {
	int x, y;
	ll t;
	Pos(int x_ = 0, int y_ = 0, ll t_ = 0) : x(x_), y(y_), t(t_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pos operator - () const { return { -x, -y }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	int Che() const { return std::max(std::abs(x), std::abs(y)); }
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
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
bool get_intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2, Pos& t, const int& F = STRONG) {
	if (!ccw(p1, p2, q1, q2)) return 0;
	if (p1.x == p2.x) { t.x = p1.x; t.y = q1.y; }
	else { t.x = q1.x; t.y = p1.y; }
	if (F == STRONG) return on_seg_strong(p1, p2, t) && on_seg_strong(q1, q2, t);
	return on_seg_strong(q1, q2, t);
}
struct Seg {
	Pos s, e;
	Seg(Pos s_ = Pos(0, 0), Pos e_ = Pos(0, 0)) : s(s_), e(e_) {}
	bool operator < (const Seg& q) const { return s.x == e.x ? s.x < q.s.x : s.y < q.s.y; }
	friend std::ostream& operator << (std::ostream& os, const Seg& S) { os << S.s << " " << S.e; return os; }
};
typedef std::vector<Seg> Segs;
Segs SV[LEN], SH[LEN];
int check(const Pos& p, const Pos& u, const Pos& d, const Seg& b, Pos& q) {
	int up = ccw(p, p + d, p + u);
	int tq1 = ccw(p, p + d, b.s);
	int tq2 = ccw(p, p + d, b.e);
	int f = tq1 * tq2;
	if (f > 0) return -1;
	Pos x1, x2;
	get_intersection(p, p + d, b.s, b.s + u, x1);
	get_intersection(p, p + d, b.e, b.e + u, x2);
	if (dot(p, p + d, x1, x2) < 0) std::swap(x1, x2);
	ll fc1 = dot(p - d, p, x1);
	ll fc2 = dot(p - d, p, x2);
	if (f < 0 || (f == 0 && (tq1 == up || tq2 == up))) {//block
		if (fc2 <= 0) return -1;
		q = x1;
		return 1;
	}
	q = Pos(fc1, fc2);
	return 0;
}
bool move(const Segs& B, Pos& p, Pos& u, Pos& d, const Pos& s) {
	static Polygon V, W;
	V.clear(); W.clear();
	Pos q;
	for (const Seg& b : B) {
		int chk = check(p, u, d, b, q);
		if (chk == -1) continue;
		else if (chk == 0) V.push_back(q);
		else W.push_back(q);
	}
	int t = 0;
	bool f = 0;
	std::sort(V.begin(), V.end());
	for (const Pos& v : V) {
		if (v.y <= 0) continue;
		if (!f) { t = v.y; f = 1; continue; }
		if (t < v.x) break;
		t = std::max(t, v.y);
	}
	Pos e = p + d * t;
	bool blk = 0;
	for (const Pos& w : W) {
		if (on_seg_strong(p, e, w)) { e = w; blk = 1; }
	}
	f = 0;
	if (on_seg_weak(p, e, s)) { f = 1; e = s; }
	int l = (p - e).Che();
	e.t = p.t + l;
	p = e;
	u *= -1;
	std::swap(u, d);
	if (blk) u *= -1, d *= -1;
	return f;
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
Polyhedron P[LEN];
Pos3D P0[LEN], NM[LEN];
ll period[LEN];
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).Euc()) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).Euc()) && dot(d1, d3, d2) > 0; }
ll dist(const Pos3D& p, const Pos3D& q) { return std::max({ std::abs(p.x - q.x), std::abs(p.y - q.y), std::abs(p.z - q.z) }); }
struct Cube {
	Pos3D a, b;
	Cube(Pos3D a_ = Pos3D(), Pos3D b_ = Pos3D()) : a(a_), b(b_) {}
	Cube& operator *= (const ll& n) { a *= n; b *= n; return *this; }
	void init() { std::cin >> *this; *this *= 2; }
	friend std::istream& operator >> (std::istream& is, Cube& c) { is >> c.a >> c.b; return is; }
	friend std::ostream& operator << (std::ostream& os, const Cube& c) { os << c.a << " " << c.b; return os; }
} C[LEN], SC[6][LEN];
typedef std::vector<Cube> Cubes;
bool cmp_x(const Cube& p, const Cube& q) { return p.a.x == q.a.x ? p.b.x < q.b.x : p.a.x < q.a.x; }
bool cmp_x_r(const Cube& p, const Cube& q) { return p.b.x == q.b.x ? p.a.x > q.a.x : p.b.x > q.b.x; }
bool cmp_y(const Cube& p, const Cube& q) { return p.a.y == q.a.y ? p.b.y < q.b.y : p.a.y < q.a.y; }
bool cmp_y_r(const Cube& p, const Cube& q) { return p.b.y == q.b.y ? p.a.y > q.a.y : p.b.y > q.b.y; }
bool cmp_z(const Cube& p, const Cube& q) { return p.a.z == q.a.z ? p.b.z < q.b.z : p.a.z < q.a.z; }
bool cmp_z_r(const Cube& p, const Cube& q) { return p.b.z == q.b.z ? p.a.z > q.a.z : p.b.z > q.b.z; }
struct Robot {
	Pos3D p, n, v, tq;
	Robot(Pos3D p_ = Pos3D(), Pos3D n_ = Pos3D(), Pos3D v_ = Pos3D()) : p(p_), n(n_), v(v_) {}
	void init() {
		std::cin >> p; p = p * 2 + _1;
		for (Pos3D* q : { &n, &v }) {
			std::cin >> dir;
			if (dir[0] == 'x') *q = Pos3D(dir[1] == '+' ? 1 : -1);
            else if (dir[0] == 'y') *q = Pos3D(0, dir[1] == '+' ? 1 : -1);
            else *q = Pos3D(0, 0, dir[1] == '+' ? 1 : -1);
		}
		p += n;
		tq = n / v;
		return;
	}
	void get_start_pos(Pos& s, Pos& u, Pos& d) const {
		if (tq.x) { s = Pos(p.y, p.z); u = Pos(n.y, n.z); d = Pos(v.y, v.z); }
		if (tq.y) { s = Pos(p.z, p.x); u = Pos(n.z, n.x); d = Pos(v.z, v.x); }
		if (tq.z) { s = Pos(p.x, p.y); u = Pos(n.x, n.y); d = Pos(v.x, v.y); }
	}
	bool slice(const Cube& c, Seg& b) const {
		if (tq.x && c.a.x <= p.x && p.x <= c.b.x) { b.s = Pos(c.a.y, c.a.z), b.e = Pos(c.b.y, c.b.z); return 1; }
		if (tq.y && c.a.y <= p.y && p.y <= c.b.y) { b.s = Pos(c.a.z, c.a.x), b.e = Pos(c.b.z, c.b.x); return 1; }
		if (tq.z && c.a.z <= p.z && p.z <= c.b.z) { b.s = Pos(c.a.x, c.a.y), b.e = Pos(c.b.x, c.b.y); return 1; }
		return 0;
	}
} R[LEN];
Pos3D convert_3D(const Pos3D& r, const Pos3D& n, const Pos& p) {
	if (n.x) return Pos3D(r.x, p.x, p.y, p.t);
	if (n.y) return Pos3D(p.y, r.y, p.x, p.t);
	return Pos3D(p.x, p.y, r.z, p.t);
}
int get_dir(const Pos3D& nm, const int& dir) {
	if (nm.x) return dir == 0 ? 2 : dir == 1 ? 3 : dir == 2 ? 4 : 5;
	if (nm.y) return dir == 0 ? 4 : dir == 1 ? 5 : dir == 2 ? 0 : 1;
	return dir;
}
int get_dir(const Pos& d) { return d.x > 0 ? 0 : d.x < 0 ? 1 : d.y > 0 ? 2 : 3; }
void get_path(const int& k) {
	P[k].reserve(10000);
	SV[k].reserve(5000);
	SH[k].reserve(5000);
	const Robot& r0 = R[k];
	P0[k] = r0.p;
	NM[k] = r0.n / r0.v;
	//if (NM[k].x) NM[k].x = P0[k].x;
	//else if (NM[k].y) NM[k].y = P0[k].y;
	//else if (NM[k].z) NM[k].z = P0[k].z;
	Pos3D n = NM[k];
	Pos s, u, d;
	r0.get_start_pos(s, u, d);
	Segs B; Seg b;
	Segs BB[4];
	for (int i = 0; i < N; i++) {
		if (r0.slice(C[i], b)) B.push_back(b);
		for (int j = 0; j < 4; j++) if (r0.slice(SC[get_dir(NM[k], j)][i], b)) BB[j].push_back(b);
	}
	//assert(B.size());
	Pos cur = s;
	Polygon L = { cur };
	while (1) {
		int dir = get_dir(d);
		const Segs& tmp = BB[dir];
		//if (move(B, cur, u, d, s)) break;
		if (move(BB[dir], cur, u, d, s)) break;
		L.push_back(cur);
	}
	int sz = L.size();
	for (int i = 0, j; i < sz; i++) {
		j = (i + 1) % sz;
		Pos3D p = convert_3D(r0.p, n, L[i]);
		P[k].push_back(p);
		Seg s = Seg(L[i], L[j]);
		bool v = L[i].x == L[j].x;
		if (v) SV[k].push_back(s);
		else SH[k].push_back(s);
	}
	std::sort(SV[k].begin(), SV[k].end());
	std::sort(SH[k].begin(), SH[k].end());
	//assert(!(cur.t % 2));
	period[k] = cur.t >> 1;
	return;
}
ll get_intersection_line(const int& i, const int& j, Pos3D& p0, Pos3D& v) {
	v = NM[i] / NM[j];
	if (!(v.x || v.y || v.z)) {
		Pos3D q = P0[i] + NM[i];
		ll fc = dot(q, P0[i], P0[j]);
		if (fc) return 0;
		int sz = P[i].size();
		for (int k = 0, l; k < sz; k++) {
			l = (k + 1) % sz;
			if (on_seg_strong(P[i][k], P[i][l], R[j].p)) {
				if (NM[i] == NM[j]) return 0;
				ll dst = dist(P[i][k], R[j].p);
				ll t = P[i][k].t + dst;
				//assert(!(t % 2));
				t >>= 1;
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
Seg projection(const Pos3D& p, const Pos3D& q, const Pos3D& n) {
	if (n.x) return Seg(Pos(p.y, p.z), Pos(q.y, q.z));
	if (n.y) return Seg(Pos(p.z, p.x), Pos(q.z, q.x));
	return Seg(Pos(p.x, p.y), Pos(q.x, q.y));
}
int intersection(const Seg& ln, const Seg& sg, const Pos3D& p0, const Pos3D& n, Pos3D& q) {
	Pos x;
	if (!get_intersection(ln.s, ln.e, sg.s, sg.e, x, WEAK)) return -1;
	int t = std::max(std::abs(sg.s.x - x.x), std::abs(sg.s.y - x.y));
	if (n.x) q = Pos3D(p0.x, x.x, x.y);
	if (n.y) q = Pos3D(x.y, p0.y, x.x);
	if (n.z) q = Pos3D(x.x, x.y, p0.z);
	t += sg.s.t;
	//assert(!(t % 2));
	return t >> 1;
}
ll collision_time(const int& k, const int& l) {
	if (R[k].p == R[l].p) return 0;
	Pos3D p0, v;
	ll f1 = get_intersection_line(k, l, p0, v);
	if (!f1) return -1;
	if (f1 < 0) return (-f1) >> 1;
	static Polygon VK, VL;
	VK.clear(); VL.clear();
	Pos3D q;
	Seg ln, sg;
	ln = projection(p0, p0 + v, NM[k]);
	const Segs& SK = ln.s.x - ln.e.x ? SV[k] : SH[k];
	int szk = SK.size();
	for (int i = 0; i < szk; i++) {
		sg = SK[i];
		int t = intersection(ln, sg, p0, NM[k], q);
		if (t >= 0) VK.push_back(Pos(v * q, t));
	}
	ln = projection(p0, p0 + v, NM[l]);
	const Segs& SL = ln.s.x - ln.e.x ? SV[l] : SH[l];
	int szl = SL.size();
	for (int i = 0; i < szl; i++) {
		sg = SL[i];
		int t = intersection(ln, sg, p0, NM[l], q);
		if (t >= 0) VL.push_back(Pos(v * q, t));	
	}
	ll ans = INF;
	int ptk = 0, ptl = 0;
    ll pk = period[k], pl = period[l];
    ExtGcd ret = ext_gcd(pk, pl);
	while (ptk < VK.size() && ptl < VL.size()) {
		if (VK[ptk].x == VL[ptl].x) {
			int h = VK[ptk].x;
			ll tk = VK[ptk].y, tl = VL[ptl].y;
			if (tk > ans) { ptk++; ptl++; continue; }
			ll diff = tl - tk;
			if (diff % ret.g == 0) {
				ll mod = pl / ret.g;
				ll n = (ret.x % mod);
				if (n < 0) n += mod;
				n = (n * ((diff / ret.g) % mod + mod) % mod) % mod;
				ll t_col = tk + n * pk;
				ans = std::min(ans, t_col);
			}
			ptk++; ptl++;
		}
		else if (VK[ptk].x < VL[ptl].x) ptk++;
		else ptl++;
	}
	return ans >= INF ? -1 : ans;
}
void cube_init() {
	std::sort(C, C + N, cmp_x); for (int i = 0; i < N; i++) SC[0][i] = C[i];
	std::sort(C, C + N, cmp_x_r); for (int i = 0; i < N; i++) SC[1][i] = C[i];
	std::sort(C, C + N, cmp_y); for (int i = 0; i < N; i++) SC[2][i] = C[i];
	std::sort(C, C + N, cmp_y_r); for (int i = 0; i < N; i++) SC[3][i] = C[i];
	std::sort(C, C + N, cmp_z); for (int i = 0; i < N; i++) SC[4][i] = C[i];
	std::sort(C, C + N, cmp_z_r); for (int i = 0; i < N; i++) SC[5][i] = C[i];
	return;
};
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(6);
	std::cin >> N >> K; ans = -1;
	for (int i = 0; i < N; i++) C[i].init();
	for (int k = 0; k < K; k++) R[k].init();
	cube_init();
	for (int k = 0; k < K; k++) get_path(k);
	for (int k = 0; k < K; k++)
		for (int l = k + 1; l < K; l++)
			if (ll t = collision_time(k, l); t >= 0)
				if (ans == -1 || t < ans) ans = t;
	if (ans < 0) std::cout << "ok\n";
	else std::cout << ans << "\n";
	return;
}
int main() { solve(); return 0; }//boj23202

/*
#include <chrono>
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(6);
	std::cin >> N >> K; ans = -1;
	for (int i = 0; i < N; i++) C[i].init();
	for (int k = 0; k < K; k++) R[k].init();

	// --- 1. 궤적 생성 (get_path) 시간 측정 시작 ---
	auto start_path = std::chrono::high_resolution_clock::now();

	for (int k = 0; k < K; k++) get_path(k);

	auto end_path = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> path_ms = end_path - start_path;

	// ---------------------------------------------
	// --- 2. 충돌 판정 (collision_time) 시간 측정 시작 ---
	auto start_collision = std::chrono::high_resolution_clock::now();

	for (int k = 0; k < K; k++)
		for (int l = k + 1; l < K; l++)
			if (ll t = collision_time(k, l); t >= 0)
				if (ans == -1 || t < ans) ans = t;

	auto end_collision = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> collision_ms = end_collision - start_collision;
	// --------------------------------------------------
	// 정답 출력 (기존 로직 유지, std::cout 사용)

	if (ans < 0) std::cout << "ok\n";
	else std::cout << ans << "\n";

	// 로그 출력 (정답 채점에 영향을 주지 않도록 std::cerr 사용)
	std::cerr << "\n========================================\n";
	std::cerr << "[실행 시간 측정 결과]\n";
	std::cerr << "- 궤적 생성 (get_path) 총 소요 시간: " << path_ms.count() << " ms\n";
	std::cerr << "- 충돌 판정 (collision) 총 소요 시간: " << collision_ms.count() << " ms\n";
	std::cerr << "========================================\n";
	return;
}
int main() { solve(); return 0; }//boj23202
*/