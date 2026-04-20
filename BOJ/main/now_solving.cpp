#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <cassert>
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
struct Pos {
	int x, y, t;
	Pos(int x_ = 0, int y_ = 0, int t_ = 0) : x(x_), y(y_), t(t_) {}
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
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
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
	bool operator < (const Seg& S) const {
		if (dot(s, e, S.s) >= 0) return 1;
		if (dot(S.s, S.e, s) >= 0) return 0;
		if (s == S.s) return (s - e).Euc() < (s - S.e).Euc();
		return (s - e).Euc() >= (S.s - e).Euc();
	}
	friend std::ostream& operator << (std::ostream& os, const Seg& S) { os << S.s << " " << S.e; return os; }
};
typedef std::vector<Seg> Segs;
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
		if (tq.x) { s = Pos(p.y, p.z); u = Pos(n.y, n.z); d = Pos(v.y, v.z); }
		if (tq.y) { s = Pos(p.z, p.x); u = Pos(n.z, n.x); d = Pos(v.z, v.x); }
		if (tq.z) { s = Pos(p.x, p.y); u = Pos(n.x, n.y); d = Pos(v.x, v.y); }
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
		if (move(B, cur, u, d, s)) break;
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
				if (NM[i] == NM[j]) return 0;
				ll dst = dist(P[i][k], R[j].p);
				ll t = P[i][k].t + dst;
				assert(!(t % 2));
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
int intersection(const Pos3D& p1, const Pos3D& p2, const Pos3D& q1, const Pos3D& q2, const Pos3D& n, Pos3D& q) {
	Pos a1, a2, b1, b2, x;
	if (n.x) { a1 = Pos(p1.y, p1.z); a2 = Pos(p2.y, p2.z); b1 = Pos(q1.y, q1.z); b2 = Pos(q2.y, q2.z); }
	if (n.y) { a1 = Pos(p1.z, p1.x); a2 = Pos(p2.z, p2.x); b1 = Pos(q1.z, q1.x); b2 = Pos(q2.z, q2.x); }
	if (n.z) { a1 = Pos(p1.x, p1.y); a2 = Pos(p2.x, p2.y); b1 = Pos(q1.x, q1.y); b2 = Pos(q2.x, q2.y); }
	if (!get_intersection(a1, a2, b1, b2, x, WEAK)) return -1;
	int t = std::max(std::abs(b1.x - x.x), std::abs(b1.y - x.y));
	if (n.x) q = Pos3D(q1.x, x.x, x.y);
	if (n.y) q = Pos3D(x.y, q1.y, x.x);
	if (n.z) q = Pos3D(x.x, x.y, q1.z);
	t += q1.t;
	assert(!(t % 2));
	return t >> 1;
}
ll collision_time(const int& k, const int& l) {
	if (R[k].p == R[l].p) return 0;
	Pos3D p0, v;
	ll f1 = get_intersection_line(k, l, p0, v);
	if (!f1) return -1;
	if (f1 < 0) return (-f1) >> 1;
	std::vector<std::pair<ull, int>> VK, VL;
	Pos3D q;
	int szk = P[k].size();
	for (int i = 0, j; i < szk; i++) {
		j = (i + 1) % szk;
		int t = intersection(p0, p0 + v, P[k][i], P[k][j], NM[k], q);
		if (t >= 0) VK.push_back({ hash(q), t });
	}
	int szl = P[l].size();
	for (int i = 0, j; i < szl; i++) {
		j = (i + 1) % szl;
		int t = intersection(p0, p0 + v, P[l][i], P[l][j], NM[l], q);
		if (t >= 0) VL.push_back({ hash(q), t });
	}
	std::sort(VK.begin(), VK.end());
	std::sort(VL.begin(), VL.end());
	ll ans = INF;
	int ptk = 0, ptl = 0;
    ll pk = period[k], pl = period[l];
    ExtGcd ret = ext_gcd(pk, pl);
	while (ptk < VK.size() && ptl < VL.size()) {
		if (VK[ptk].first == VL[ptl].first) {
			ull h = VK[ptk].first;
			ll tk = VK[ptk].second;
			ll tl = VL[ptl].second;
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
		else if (VK[ptk].first < VL[ptl].first) ptk++;
		else ptl++;
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
	for (int k = 0; k < K; k++) {
		for (int l = k + 1; l < K; l++) {
			ll t = collision_time(k, l);
			if (t >= 0) {
				if (ans == -1 || t < ans) ans = t;
			}
		}
	}
	if (ans < 0) std::cout << "ok\n";
	else std::cout << ans << "\n";
	return;
}
int main() { solve(); return 0; }//boj23202
