#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <random>
#include <array>
#include <tuple>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 1e5 + 1;
const ld TOL = 1e-7;
const ll MOD = 1e9 + 7;
const ld PI = acos(-1);
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo, const ld& hi) { return std::min(hi, std::max(lo, x)); }
ld flip(ld lat) {
	if (zero(lat - PI * .5) || zero(lat + PI * .5)) return 0;
	if (zero(lat)) return PI * .5;
	if (lat > 0) return PI * .5 - lat;
	if (lat < 0) return -(PI * .5) - lat;
	return INF;
}
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
//ll pow_fuck(ll a, ll b) {
//	ll ret = 1;
//	while (b) {
//		if (b & 1) ret = ret * a % MOD;
//		a = a * a % MOD;
//		b >>= 1;
//	}
//	return ret;
//}
//ll powmod(ll a, ll b) {
//	ll res = 1; a %= MOD;
//	assert(b >= 0);
//	for (; b; b >>= 1) {
//		if (b & 1) res = res * a % MOD;
//		a = a * a % MOD;
//	}
//	return res;
//}
struct Info { ll area, l, r; };

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

//freopen("../../../input_data/triathlon_tests/triath.20", "r", stdin);
//freopen("../../../input_data/triathlon_tests/triathlon_out.txt", "w", stdout);

//Euler characteristic : v - e + f == 1
//Pick`s Theorem : A = i + b/2 - 1

//2D============================================================================//

int N, M, T, Q;
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
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
//bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
bool cmpt(const Pos& p, const Pos& q) {
	bool f0 = O < p;
	bool f1 = O < q;
	if (f0 != f1) return f0;
	ll tq = p / q;
	return !tq ? p.Euc() < q.Euc() : tq > 0;
}
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
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2, const int& f = STRONG) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	if (f == WEAK) return f1 && f2;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
void get_area_memo(Pos H[], ll memo[], const int& sz) {
	memo[0] = 0;
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		memo[i + 1] = cross(Pos(0, 0), cur, nxt) + memo[i];//memo[sz] == convex hull's area
	}
	return;
}
void get_round_memo(Polygon& H, ld memo[]) {
	int sz = H.size();
	memo[0] = .0;
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		memo[i + 1] = (cur - nxt).mag() + memo[i];//memo[sz] == convex hull's round
	}
	return;
}
ll area(Pos H[], const int& sz) {
	ll a = 0;
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
ll area(Polygon& H) {
	ll a = 0; int sz = H.size();
	for (int i = 0; i < sz; i++) a += H[i] / H[(i + 1) % sz];
	return a;
}
void norm(Polygon& H) { if (area(H) < 0) std::reverse(H.begin(), H.end()); }
bool inner_check(Pos p0, Pos p1, Pos p2, const Pos& t) {
	if (ccw(p0, p1, p2) < 0) std::swap(p1, p2);
	return ccw(p0, p1, t) >= 0 && ccw(p1, p2, t) >= 0 && ccw(p2, p0, t) >= 0;
}
int inner_check(Pos H[], const int& sz, const Pos& p) {//concave
	int cnt{ 0 };
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (cur.y == nxt.y) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y <= p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return (cnt & 1) * 2;
}
int inner_check(const Polygon& H, const Pos& p) {//concave
	int cnt = 0, sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos cur = H[i], nxt = H[(i + 1) % sz];
		if (on_seg_strong(cur, nxt, p)) return 1;
		if (cur.y == nxt.y) continue;
		if (nxt.y < cur.y) std::swap(cur, nxt);
		if (nxt.y <= p.y || cur.y > p.y) continue;
		cnt += ccw(cur, nxt, p) > 0;
	}
	return (cnt & 1) * 2;
}
int inner_check_bi_search(Pos H[], const int& sz, const Pos& p) {//convex
	if (!sz) return -1;
	if (sz == 1) return p == H[0] ? 0 : -1;
	if (sz == 2) return on_seg_strong(H[0], H[1], p) ? 0 : -1;
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return -1;
	if (on_seg_strong(H[0], H[1], p) || on_seg_strong(H[0], H[sz - 1], p)) return 0;
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	if (cross(H[s], H[e], p) > 0) return 1;
	else if (on_seg_strong(H[s], H[e], p)) return 0;
	else return -1;
}
int inner_check_bi_search(const Polygon& H, const Pos& p) {//convex
	int sz = H.size();
	if (!sz) return -1;
	if (sz == 1) return p == H[0] ? 0 : -1;
	if (sz == 2) return on_seg_strong(H[0], H[1], p) ? 0 : -1;
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return -1;
	if (on_seg_strong(H[0], H[1], p) || on_seg_strong(H[0], H[sz - 1], p)) return 0;
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	if (cross(H[s], H[e], p) > 0) return 1;
	else if (on_seg_strong(H[s], H[e], p)) return 0;
	else return -1;
}
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	std::sort(C.begin(), C.end());
	C.erase(unique(C.begin(), C.end()), C.end());
	if (C.size() <= 2) { for (const Pos& pos : C) H.push_back(pos); }
	else {
		for (int i = 0; i < C.size(); i++) {
			while (H.size() > 1 && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
		int s = H.size() + 1;
		for (int i = C.size() - 1; i >= 0; i--) {
			while (H.size() > s && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) <= 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
	}
	return H;
}
Polygon graham_scan(Polygon& C) {
	Polygon H;
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
ll rotating_calipers(const Polygon& H) {
	int sz = H.size();
	if (sz < 3) return -1;
	auto jaw = [&](const int& i, const int& f) -> ll {
		return (H[(i + 1) % sz] - H[i]) / (H[(f + 1) % sz] - H[f]);
		};
	ll ret = 0;
	for (int i = 0, j = 1; i < sz; i++) {
		while (jaw(i, j) > 0) {
			ret = std::max(ret, (H[i] - H[j]).Euc());
			j = (j + 1) % sz;
		}
		ret = std::max(ret, (H[i] - H[j]).Euc());
	}
	return ret;
}
Info find_tangent_bi_search(Pos H[], const int& sz, const Pos& p) {
	int i1{ 0 }, i2{ 0 };
	int ccw1 = ccw(p, H[0], H[1]), ccwN = ccw(p, H[0], H[sz - 1]);
	if (ccw1 * ccwN >= 0) {
		i1 = 0;
		if (!ccw1 && dot(p, H[1], H[0]) > 0) i1 = 1;
		if (!ccwN && dot(p, H[sz - 1], H[0]) > 0) i1 = sz - 1;
		int s = 0, e = sz - 1, m;
		if (!ccw1) s += 1;
		if (!ccwN) e -= 1;
		bool f = ccw(p, H[s], H[s + 1]) >= 0;
		while (s < e) {
			m = s + e >> 1;
			Pos p1 = p, cur = H[m], nxt = H[(m + 1) % sz];
			if (!f) std::swap(p1, cur);//normalize
			if (ccw(p1, cur, nxt) > 0) s = m + 1;
			else e = m;
		}
		i2 = s;
		if (!ccw(p, H[i2], H[(i2 + 1) % sz]) && dot(p, H[(i2 + 1) % sz], H[i2]) > 0) i2 = (i2 + 1) % sz;
	}
	else {
		//divide hull
		int s = 0, e = sz - 1, k, m;
		bool f = ccw1 > 0 && ccwN < 0;//if H[k] is between H[0] && p
		while (s + 1 < e) {
			k = s + e >> 1;
			int CCW = ccw(H[0], H[k], p);
			if (!f) CCW *= -1;//normailze
			if (CCW > 0) s = k;
			else e = k;
		}

		//search lower hull
		int s1 = 0, e1 = s;
		while (s1 < e1) {
			m = s1 + e1 >> 1;
			Pos p1 = p, cur = H[m], nxt = H[(m + 1) % sz];
			if (!f) std::swap(p1, cur);//normalize
			if (ccw(p1, cur, nxt) > 0) s1 = m + 1;
			else e1 = m;
		}
		i1 = s1;
		if (!ccw(p, H[i1], H[(i1 + 1) % sz]) && dot(p, H[(i1 + 1) % sz], H[i1]) > 0) i1 = (i1 + 1) % sz;

		//search upper hull
		int s2 = e, e2 = sz - 1;
		while (s2 < e2) {
			m = s2 + e2 >> 1;
			Pos p1 = p, cur = H[m], nxt = H[(m + 1) % sz];
			if (!f) std::swap(p1, cur);//normalize
			if (ccw(p1, cur, nxt) < 0) s2 = m + 1;
			else e2 = m;
		}
		i2 = s2;
		if (!ccw(p, H[i2], H[(i2 + 1) % sz]) && dot(p, H[(i2 + 1) % sz], H[i2]) > 0) i2 = (i2 + 1) % sz;
	}
	if (i2 < i1) std::swap(i2, i1);//normalize
	return { 0, i2, i1 };
}
Info get_inner_area(Pos H[], ll memo[], const int& sz, const Pos& p) {
	Info tangent = find_tangent_bi_search(H, sz, p);
	ll i1 = tangent.r, i2 = tangent.l;
	ll tri = cross(O, H[i1], H[i2]);
	ll area = memo[i2] - memo[i1] - tri;
	if (cross(p, H[i1], H[i2]) < 0) area = memo[sz] - area, std::swap(i1, i2);//normalize
	area += std::abs(cross(p, H[i1], H[i2]));
	return { area, i2, i1 };
}
ld find_inx_get_area_bi_search(Pos H_in[], ll memo_in[], const int& sz_in, Pos H_out[], ll memo_out[], const int& sz_out, const Pos& p) {
	Info info = get_inner_area(H_in, memo_in, sz_in, p);
	Pos vr = H_in[info.r], vl = H_in[info.l];
	int ir, il;
	ld wing_r{ 0 }, wing_l{ 0 };

	//divide hull
	int s = 0, e = sz_out - 1, k, m;
	while (s + 1 < e) {
		k = s + e >> 1;
		int CCW = ccw(H_out[0], H_out[k], p);
		if (CCW > 0) s = k;
		else e = k;
	}
	Pos S = H_out[s], E = H_out[e];

	//find r-intersection
	int sr{ 0 }, er{ 0 };
	Pos SR, ER;
	if (ccw(p, S, vr) >= 0 && ccw(p, E, vr) <= 0) sr = s, er = e;//if vr is in p-S-E tri.
	else {
		if (ccw(H_out[0], p, vr) > 0) sr = e, er = sz_out;
		if (ccw(H_out[0], p, vr) < 0) sr = 0, er = s;
		while (sr + 1 < er) {
			m = sr + er >> 1;
			int CCW = ccw(p, H_out[m % sz_out], vr);
			if (CCW > 0) sr = m;
			else er = m;
		}
	}
	SR = H_out[sr % sz_out], ER = H_out[er % sz_out];
	ir = er % sz_out;
	ll trir = std::abs(cross(p, SR, ER));
	ll ar = std::abs(cross(p, vr, SR));
	ll br = std::abs(cross(p, vr, ER));
	wing_r = trir * (ld)br / (ar + br);
	if (!cross(p, vr, H_out[er % sz_out])) wing_r = 0;

	//find l-intersection
	int sl{ 0 }, el{ 0 };
	Pos SL, EL;
	if (ccw(p, S, vl) >= 0 && ccw(p, E, vl) <= 0) sl = s, el = e;//if vl is in p-S-E tri.
	else {
		if (ccw(H_out[0], p, vl) > 0) sl = e, el = sz_out;
		if (ccw(H_out[0], p, vl) < 0) sl = 0, el = s;
		while (sl + 1 < el) {
			m = sl + el >> 1;
			int CCW = ccw(p, H_out[m % sz_out], vl);
			if (CCW > 0) sl = m;
			else el = m;
		}
	}
	SL = H_out[sl % sz_out], EL = H_out[el % sz_out];
	il = sl % sz_out;
	ll tril = std::abs(cross(p, SL, EL));
	ll al = std::abs(cross(p, vl, SL));
	ll bl = std::abs(cross(p, vl, EL));
	wing_l = tril * (ld)al / (al + bl);
	if (!cross(p, vl, H_out[sl % sz_out])) wing_l = 0;

	//DEBUG
	//std::cout << "in R: " << info.r << " in L: " << info.l << " out R: " << ir << " out L: " << il << "\n";
	//std::cout << "wing R: " << wing_r << " wing L: " << wing_l << "\n";
	//std::cout << "wing R: " << trir * (ld)ar / (ar + br) << " wing L: " << tril * (ld)bl / (al + bl) << "\n";
	//std::cout << "inner: " << info.area << "\n";

	//get_shadow
	ld area{ 0 };
	if (sr == sl) {//if 2 intersections on the same segment
		area = -(ld)(info.area + std::abs(cross(p, H_out[ir], H_out[il]))) + (wing_r + wing_l);
	}
	else {
		bool f = ir > il;
		if (ir > il) std::swap(ir, il);//normalize
		ll tri = cross(O, H_out[ir], H_out[il]);
		ll tmp = memo_out[il] - memo_out[ir] - tri;
		if (f) tmp = memo_out[sz_out] - tmp;
		area = -(ld)(info.area - tmp - std::abs(cross(p, H_out[ir], H_out[il]))) + (wing_r + wing_l);
	}
	return area * .5;
}
Polygon minkowski_sum(Polygon& A, Polygon& B, const bool& f = 0) {
	int sza = A.size();
	int szb = B.size();
	Polygon MS;

	if (f) { for (Pos& p : A) p *= 2; for (Pos& p : B) p *= -1; }
	int a = 0, b = 0;
	for (int i = 0; i < sza; i++) if (A[i] < A[a]) a = i;
	for (int j = 0; j < szb; j++) if (B[j] < B[b]) b = j;
	MS.push_back(A[a] + B[b]);
	//std::cout << sza + szb << "\n";
	for (int k = 1; k < sza + szb; k++) {
		//std::cout << "DEBUG\n";
		Pos a1 = A[(a + 1) % sza] + B[b];
		Pos b1 = B[(b + 1) % szb] + A[a];
		if (ccw(MS.back(), a1, b1) > 0) MS.push_back(a1), a = (a + 1) % sza;// , std::cout << a1;
		else MS.push_back(b1), b = (b + 1) % szb;// , std::cout << b1;
	}

	if (f) { for (Pos& p : A) p /= 2; for (Pos& p : B) p *= -1; }
	return graham_scan(MS);
}
Pos inner_check_bi_search_seg(const Polygon& H, const Pos& p) {//convex
	int sz = H.size();
	if (!sz) return Pos(-1, -1);
	if (sz == 1) return p == H[0] ? Pos(0, 0) : Pos(-1, -1);
	if (sz == 2) {
		int i1 = -1, i2 = -1;
		if (H[0] == p) i1 = 0;
		if (H[1] == p) i2 = 1;
		return Pos(i1, i2);
	}
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return Pos(-1, -1);
	if (H[0] == p) return Pos(sz - 1, 1);
	if (H[1] == p) return Pos(0, 2 % sz);
	if (H[sz - 1] == p) return Pos(sz - 2, 0);
	if (on_seg_weak(H[0], H[1], p)) return Pos(0, 1);
	if (on_seg_weak(H[0], H[sz - 1], p)) return Pos(sz - 1, 0);
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	//std::cout << "DEBUG:: " << s << " " << e << "\n";
	if (cross(H[s], H[e], p) < 0) return Pos(-1, -1);
	if (H[s] == p) return Pos((s - 1 + sz) % sz, e);
	if (H[e] == p) return Pos(s, (e + 1) % sz);
	if (on_seg_weak(H[s], H[e], p)) return Pos(s, e);
	//std::cout << "DEBUG:: fuck\n";
	return Pos(sz + 1, sz + 1);
}
Pos find_tangent_bi_search(const Polygon& H, const Pos& p) {
	int sz = H.size();
	Pos IN = Pos(sz + 1, sz + 1);
	Pos F = inner_check_bi_search_seg(H, p);
	//std::cout << "inner_bi:: " << F << "\n";
	if (F == IN) return INVAL;
	if (F != INVAL) return F;
	int i1{ 0 }, i2{ 0 };
	int ccw1 = ccw(p, H[0], H[1]), ccwN = ccw(p, H[0], H[sz - 1]);
	if (ccw1 * ccwN >= 0) {
		i1 = 0;
		if (!ccw1 && dot(p, H[1], H[0]) > 0) i1 = 1;
		if (!ccwN && dot(p, H[sz - 1], H[0]) > 0) i1 = sz - 1;
		int s = 0, e = sz - 1, m;
		if (!ccw1) s += 1;
		if (!ccwN) e -= 1;
		bool f = ccw(p, H[s], H[s + 1]) >= 0;
		while (s < e) {
			m = s + e >> 1;
			const Pos& p1 = p, & cur = H[m], & nxt = H[(m + 1) % sz];
			int CCW = ccw(p1, cur, nxt);
			if (!f) CCW *= -1;//normailze
			if (CCW > 0) s = m + 1;
			else e = m;
		}
		i2 = s;
		if (!ccw(p, H[i2], H[(i2 + 1) % sz]) && dot(p, H[(i2 + 1) % sz], H[i2]) > 0) i2 = (i2 + 1) % sz;
	}
	else {
		//divide hull
		int s = 0, e = sz - 1, k, m;
		bool f = ccw1 > 0 && ccwN < 0;//if H[k] is between H[0] && p
		while (s + 1 < e) {
			k = s + e >> 1;
			int CCW = ccw(H[0], H[k], p);
			if (!f) CCW *= -1;//normailze
			if (CCW > 0) s = k;
			else e = k;
		}

		//search lower hull
		int s1 = 0, e1 = s;
		while (s1 < e1) {
			m = s1 + e1 >> 1;
			const Pos& p1 = p, & cur = H[m], & nxt = H[(m + 1) % sz];
			int CCW = ccw(p1, cur, nxt);
			if (!f) CCW *= -1;//normailze
			if (CCW > 0) s1 = m + 1;
			else e1 = m;
		}
		i1 = s1;
		if (!ccw(p, H[i1], H[(i1 + 1) % sz]) && dot(p, H[(i1 + 1) % sz], H[i1]) > 0) i1 = (i1 + 1) % sz;

		//search upper hull
		int s2 = e, e2 = sz - 1;
		while (s2 < e2) {
			m = s2 + e2 >> 1;
			const Pos& p1 = p, & cur = H[m], & nxt = H[(m + 1) % sz];
			int CCW = ccw(p1, cur, nxt);
			if (!f) CCW *= -1;//normailze
			if (CCW < 0) s2 = m + 1;
			else e2 = m;
		}
		i2 = s2;
		if (!ccw(p, H[i2], H[(i2 + 1) % sz]) && dot(p, H[(i2 + 1) % sz], H[i2]) > 0) i2 = (i2 + 1) % sz;
	}
	if (ccw(p, H[i1], H[i2]) < 0) std::swap(i1, i2);
	return Pos(i2, i1);
}

struct Vec {
	ll vy, vx;
	Vec(ll y_ = 0, ll x_ = 0) : vy(y_), vx(x_) {}
	bool operator == (const Vec& v) const { return vy == v.vy && vx == v.vx; }
	bool operator < (const Vec& v) const { return vy == v.vy ? vx < v.vx : vy < v.vy; }
	ll operator / (const Vec& v) const { return vy * v.vx - vx * v.vy; }
}; const Vec Zero = { 0, 0 };

struct Line {//ax + by = c
	Vec s;
	ll c;
	Line(Vec v_ = Vec(0, 0), ll c_ = 0) : s(v_), c(c_) {}
	Line(Vec v_ = Vec(0, 0), Pos p = Pos(0, 0)) : s(v_) { c = s.vy * p.x + s.vx * p.y; }
	Line(Pos ps = Pos(0, 0), Pos pe = Pos(0, 0)) {
		ll dy, dx, c;
		dy = pe.y - ps.y;
		dx = ps.x - pe.x;
		ll _gcd = gcd(std::abs(dy), std::abs(dx));
		dy /= _gcd; dx /= _gcd;
		s = Vec(dy, dx);
		c = dy * ps.x + dx * ps.y;
	}
	bool operator == (const Line& l) const { return s == l.s && c == l.c; }
	bool operator < (const Line& l) const {
		bool f1 = Zero < s;
		bool f2 = Zero < l.s;
		if (f1 != f2) return f1;
		ll CCW = s / l.s;
		return !CCW ? c < l.c : CCW > 0;
	}
	friend std::ostream& operator << (std::ostream& os, const Line& l) { os << l.s.vy << " " << l.s.vx << " " << l.c;return os; }
};
Line L(const Pos& s, const Pos& e) {
	ll dy, dx, c;
	dy = e.y - s.y;
	dx = s.x - e.x;
	ll _gcd = gcd(std::abs(dy), std::abs(dx));
	dy /= _gcd; dx /= _gcd;
	c = dy * s.x + dx * s.y;
	return Line(Vec(dy, dx), c);
}

struct Seg {
	Line l;
	Pos s, e;
	Seg(Line l_ = Line(Vec(0, 0), 0), Pos s_ = Pos(0, 0), Pos e_ = Pos(0, 0)) : l(l_), s(s_), e(e_) {}
	bool operator == (const Seg& S) const { return l == S.l && s == S.s && e == S.e; }
	bool operator != (const Seg& S) const { return !(*this == S); }
	bool operator < (const Seg& S) const { return (l == S.l) ? (s == S.s) ? e < S.e : s < S.s : l < S.l; }
	friend std::ostream& operator << (std::ostream& os, const Seg& S) { os << "DEBUG::Seg l: " << S.l << " | s: " << S.s << " | e: " << S.e << " DEBUG::Seg\n";return os; }
};
void make_seg(std::vector<Seg>& V, const Pos& x1, const Pos& x2) {
	Pos d1 = x1, d2 = x2;
	assert(d2 != d1);
	if (d2 < d1) std::swap(d1, d2);
	V.push_back(Seg(L(d1, d2), d1, d2));
	return;
}
void sweep(std::vector<Seg>& V, std::vector<Seg>& V2) {
	std::sort(V.begin(), V.end());
	Pos tmp;
	int sz = V.size();
	for (int i = 0, j; i < sz; i = j) {//sweeping
		j = i;
		while (j < sz && V[i].l == V[j].l) j++;
		for (int k = i; k < j - 1; k++) {
			int nxt = k + 1;
			if (V[k].e < V[nxt].s) continue;
			else if (V[k].e == V[nxt].s) {
				tmp = V[k].s;
				V[k].s = V[nxt].s;
				V[nxt].s = tmp;
			}
			else if (V[nxt].e < V[k].e) {
				tmp = V[k].e;
				V[k].e = V[nxt].s;
				V[nxt].s = V[nxt].e;
				V[nxt].e = tmp;
			}
			else if (V[k].e <= V[nxt].e) {
				tmp = V[k].e;
				V[k].e = V[nxt].s;
				V[nxt].s = tmp;
			}
		}
		for (int k = i; k < j; k++) if (V[k].s != V[k].e) V2.push_back(V[k]);
	}
	//std::sort(V2.begin(), V2.end());
	return;
}

//3D============================================================================//

struct Pos3D {
	ll x, y, z;
	Pos3D(ll x_ = 0, ll y_ = 0, ll z_ = 0) : x(x_), y(y_), z(z_) {}
	bool operator == (const Pos3D& p) const { return x == p.x && y == p.y && z == p.z; }
	bool operator != (const Pos3D& p) const { return x != p.x || y != p.y || z != p.z; }
	bool operator < (const Pos3D& p) const { return x == p.x ? y == p.y ? z < p.z : y < p.y : x < p.x; }
	ll operator * (const Pos3D& p) const { return x * p.x + y * p.y + z * p.z; }
	Pos3D operator / (const Pos3D& p) const {
		Pos3D ret;
		ret.x = y * p.z - z * p.y;
		ret.y = z * p.x - x * p.z;
		ret.z = x * p.y - y * p.x;
		return ret;
	}
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	Pos3D operator * (const ll& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ll& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ll& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ll& n) { x /= n; y /= n; z /= n; return *this; }
	ll Euc() const { return x * x + y * y + z * z; }
	ld mag() const { return sqrtl(Euc()); }
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z << "\n"; return os; }
} candi[LEN];
const Pos3D O3D = { 0, 0, 0 };
const Pos3D X_axis = { 1, 0, 0 };
const Pos3D Y_axis = { 0, 1, 0 };
const Pos3D Z_axis = { 0, 0, 1 };
const Pos3D MAXP3D = { INF, INF, INF };
typedef std::vector<Pos3D> Polygon3D;
typedef std::vector<Polygon3D> Polyhedron;
using Face = std::array<int, 3>;
std::vector<Face> Hull3D;
struct Edge {
	int face_num, edge_num;
	Edge(int t = 0, int v = 0) : face_num(t), edge_num(v) {}
};
struct Line3D {
	Pos3D dir, p0;
	Line3D(Pos3D d = Pos3D(), Pos3D p = Pos3D()) : dir(d), p0(p) {}
};
struct Planar {
	Pos3D norm, p0;
	Planar(Pos3D n = Pos3D(), Pos3D p = Pos3D()) : norm(n), p0(p) {}
	friend std::istream& operator >> (std::istream& is, Planar& P) { is >> P.norm >> P.p0; return is; }
	friend std::ostream& operator << (std::ostream& os, const Planar& P) { os << P.norm << " " << P.p0; return os; }
};
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return !zero(cross(d1, d2, d3).mag()) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && dot(d1, d3, d2) > 0; }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
ld area(const std::vector<Pos3D>& H, const Pos3D& norm) {
	if (H.size() < 3) return 0;
	ll ret = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		Pos3D cur = H[i], nxt = H[(i + 1) % sz];
		ret += cross(H[0], cur, nxt) * norm / norm.Euc();
	}
	return std::abs(ret * .5);
}
//std::vector<Pos3D> graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
ld graham_scan(std::vector<Pos3D>& C, const Pos3D& norm) {
	//if (C.size() < 3) {
	//	std::sort(C.begin(), C.end());
	//	return C;
	// }
	if (C.size() < 3) return 0;
	std::vector<Pos3D> H;
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pos3D& p, const Pos3D& q) -> bool {
		int ret = ccw(C[0], p, q, norm);
		if (!ret) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i], norm) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	//return H;
	return area(H, norm);
}
bool collinear(const Pos3D& a, const Pos3D& b, const Pos3D& c) { return zero(((b - a) / (c - b)).Euc()); }
bool coplanar(const Pos3D& a, const Pos3D& b, const Pos3D& c, const Pos3D& p) { return zero(cross(a, b, c) * (p - a)); }
bool above(const Pos3D& a, const Pos3D& b, const Pos3D& c, const Pos3D& p) { return cross(a, b, c) * (p - a) > 0; }
int prep(std::vector<Pos3D>& p) {//refer to Koosaga'
	shuffle(p.begin(), p.end(), std::mt19937(0x14004));
	int dim = 1;
	for (int i = 1; i < p.size(); i++) {
		if (dim == 1) {
			if (p[0] != p[i])
				std::swap(p[1], p[i]), ++dim;
		}
		else if (dim == 2) {
			if (!collinear(p[0], p[1], p[i]))
				std::swap(p[2], p[i]), ++dim;
		}
		else if (dim == 3) {
			if (!coplanar(p[0], p[1], p[2], p[i]))
				std::swap(p[3], p[i]), ++dim;
		}
	}
	assert(dim == 4);
	return dim;
}
struct Face {
	int v[3];
	Face(int a = 0, int b = 0, int c = 0) { v[0] = a; v[1] = b; v[2] = c; }
	Pos3D norm(std::vector<Pos3D>& C) const { return cross(C[v[0]], C[v[1]], C[v[2]]); }
	Planar P(std::vector<Pos3D>& C) const { return Planar(norm(C), C[v[0]]); }
};
ld dist(const std::vector<Pos3D>& C, const Face& F, const Pos3D& p) {
	Pos3D norm = cross(C[F.v[0]], C[F.v[1]], C[F.v[2]]);
	ll ret = norm * (p - C[F.v[0]]);
	return std::abs(ret / (norm.mag()));
}
ld get_min_dist(const std::vector<Pos3D>& C, const std::vector<Face>& F, const Pos3D& p) {
	ld MIN = INF;
	for (const Face& face : F) MIN = std::min(MIN, dist(C, face, p));
	return MIN;
}
struct Edge {
	int face_num, edge_num;
	Edge(int t = 0, int v = 0) : face_num(t), edge_num(v) {}
};
bool col = 0, cop = 0;
std::vector<Face> convex_hull_3D(std::vector<Pos3D>& candi) {//incremental construction
	// 3D Convex Hull in O(n log n)
	// Very well tested. Good as long as not all points are coplanar
	// In case of collinear faces, returns arbitrary triangulation
	// Credit: Benq
	// refer to Koosaga'
	col = 0, cop = 0;
	int suf = prep(candi);
	if (suf <= 2) { col = 1; return {}; };
	if (suf == 3) { cop = 1; return {}; };
	int sz = candi.size();
	std::vector<Face> faces;
	std::vector<int> active;//whether face is active - face faces outside 
	std::vector<std::vector<int>> vis(sz);//faces visible from each point
	std::vector<std::vector<int>> rvis;//points visible from each face
	std::vector<std::array<Edge, 3>> other;//other face adjacent to each edge of face
	auto ad = [&](const int& a, const int& b, const int& c) -> void {//add face
		faces.push_back({ a, b, c });
		active.push_back(1);
		rvis.emplace_back();
		other.emplace_back();
		return;
		};
	auto visible = [&](const int& a, const int& b) -> void {
		vis[b].push_back(a);
		rvis[a].push_back(b);
		return;
		};
	auto abv = [&](const int& a, const int& b) -> bool {//above
		Face tri = faces[a];
		return above(candi[tri.v[0]], candi[tri.v[1]], candi[tri.v[2]], candi[b]);
		};
	auto edge = [&](const Edge& e) -> pi {
		return { faces[e.face_num].v[e.edge_num], faces[e.face_num].v[(e.edge_num + 1) % 3] };
		};
	auto glue = [&](const Edge& a, const Edge& b) -> void {//link two faces by an edge
		pi x = edge(a); assert(edge(b) == pi(x.second, x.first));
		other[a.face_num][a.edge_num] = b;
		other[b.face_num][b.edge_num] = a;
		return;
		};//ensure face 0 is removed when i = 3

	ad(0, 1, 2), ad(0, 2, 1);
	if (abv(1, 3)) std::swap(candi[1], candi[2]);
	for (int i = 0; i < 3; i++) glue({ 0, i }, { 1, 2 - i });
	for (int i = 3; i < N; i++) visible(abv(1, i), i);//coplanar points go in rvis[0]

	std::vector<int> label(sz, -1);
	for (int i = 3; i < sz; i++) {//incremental construction
		std::vector<int> rem;
		for (auto& v : vis[i]) if (active[v]) { active[v] = 0, rem.push_back(v); }
		if (!rem.size()) continue;//hull unchanged

		int st = -1;//start idx
		for (const int& v : rem) {
			for (int j = 0; j < 3; j++) {
				int o = other[v][j].face_num;
				if (active[o]) {//create new face!
					int idx1, idx2;
					std::tie(idx1, idx2) = edge({ v, j });
					ad(idx1, idx2, i);
					st = idx1;
					int cur = rvis.size() - 1;
					label[idx1] = cur;

					std::vector<int> tmp;
					set_union(rvis[v].begin(), rvis[v].end(), rvis[o].begin(), rvis[o].end(), back_inserter(tmp));
					//merge sorted vectors ignoring duplicates

					for (auto& x : tmp) if (abv(cur, x)) visible(cur, x);
					//if no rounding errors then guaranteed that only x > i matters

					glue({ cur, 0 }, other[v][j]);//glue old, new face
				}
			}
		}
		for (int x = st, y; ; x = y) {//glue new faces together
			int X = label[x];
			glue({ X, 1 }, { label[y = faces[X].v[1]], 2 });
			if (y == st) break;
		}
	}
	std::vector<Face> hull3D;
	for (int i = 0; i < faces.size(); i++) if (active[i]) hull3D.push_back(faces[i]);
	return hull3D;
}

//=================================================//

using pi = std::pair<int, int>;
#define sz(v) ((int)((v).size()))
const int MAXN = 1005;

struct vec3 {
	ll x, y, z;
	vec3() : x(0), y(0), z(0) {}
	vec3(ll a, ll b, ll c) : x(a), y(b), z(c) {}
	vec3 operator*(const vec3& v) const { return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	vec3 operator-(const vec3& v) const { return vec3(x - v.x, y - v.y, z - v.z); }
	vec3 operator-() const { return vec3(-x, -y, -z); }
	ll dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
};

struct twoset {
	int a, b;
	void insert(int x) { (a == -1 ? a : b) = x; }
	bool contains(int x) { return a == x || b == x; }
	void erase(int x) { (a == x ? a : b) = -1; }
	int size() { return (a != -1) + (b != -1); }
} E[MAXN][MAXN]; // i < j

struct face {
	vec3 norm;
	ll disc;
	int I[3];
};

face make_face(int i, int j, int k, int ii, std::vector<vec3>& A) {// p^T * norm < disc
	E[i][j].insert(k); E[i][k].insert(j); E[j][k].insert(i);
	face f; f.I[0] = i, f.I[1] = j, f.I[2] = k;
	f.norm = (A[j] - A[i]) * (A[k] - A[i]);
	f.disc = f.norm.dot(A[i]);
	if (f.norm.dot(A[ii]) > f.disc) {
		f.norm = -f.norm;
		f.disc = -f.disc;
	}
	return f;
}

std::vector<face> get_hull(std::vector<vec3>& A) {
	int N = A.size();
	std::vector<face> faces; memset(E, -1, sizeof(E));
	faces.push_back(make_face(0, 1, 2, 3, A));
	faces.push_back(make_face(0, 1, 3, 2, A));
	faces.push_back(make_face(0, 2, 3, 1, A));
	faces.push_back(make_face(1, 2, 3, 0, A));
	for (int i = 4; i < N; ++i) {
		for (int j = 0; j < faces.size(); ++j) {
			face f = faces[j];
			if (f.norm.dot(A[i]) > f.disc) {
				E[f.I[0]][f.I[1]].erase(f.I[2]);
				E[f.I[0]][f.I[2]].erase(f.I[1]);
				E[f.I[1]][f.I[2]].erase(f.I[0]);
				faces[j--] = faces.back();
				faces.pop_back();
			}
		}
		int nf = faces.size();
		for (int j = 0; j < nf; ++j) {
			face f = faces[j];
			for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) {
				int c = 3 - a - b;
				if (E[f.I[a]][f.I[b]].size() == 2) continue;
				faces.push_back(make_face(f.I[a], f.I[b], i, f.I[c], A));
			}
		}
	}
	return faces;
}

void solve() {
	int n; scanf("%d", &n);
	if (n == 0) exit(0);
	std::vector<vec3> inp;
	for (int i = 1; i <= n; i++) {
		int x, y, z; scanf("%d %d %d", &x, &y, &z);
		inp.emplace_back(x, y, z);
	}
	auto h = get_hull(inp);
	int q; scanf("%d", &q);
	while (q--) {
		int x, y, z; scanf("%d %d %d", &x, &y, &z);
		vec3 qry(x, y, z);
		double ret = 1e18;
		for (auto& i : h) {
			double dist = i.disc - i.norm.dot(qry);
			dist = dist * dist / i.norm.dot(i.norm);
			ret = std::min(ret, dist);
		}
		printf("%.4f\n", sqrt(ret));
	}
}

int main() {
	int tc = 1e9;// scanf("%d",&tc);
	for (int i = 1; i <= tc; i++) {
		solve();
	}
}