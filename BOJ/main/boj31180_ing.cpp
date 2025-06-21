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
#include <deque>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 5e5 + 10;
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
ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define LO x
#define HI y

#define LINE 1
#define CIRCLE 2

#define STRONG 0
#define WEAK 1

//freopen("../../../input_data/triathlon_tests/triath.20", "r", stdin);
//freopen("../../../input_data/triathlon_tests/triathlon_out.txt", "w", stdout);

int N, M, T, Q;
struct Pos {
	int x, y;
	int i, f;
	//ll x, y;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1, int f_ = -1) : x(x_), y(y_), i(i_), f(f_) {}
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
}; const Pos O = Pos(0, 0, -1, -1);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpi(const Pos& p, const Pos& q) { return p.i < q.i; }
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
Polygon graham_scan(Polygon& C, const int& f = 0) {
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
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i]) <= f)
			H.pop_back();
		H.push_back(C[i]);
	}
	return H;
}
ld R[LEN];
int RI[LEN], LI[LEN];
int RC[LEN], LC[LEN];
ld RR[LEN], LR[LEN];
ld yes_O(Polygon& P, const ld& r_) {
	std::cout << "YES::\n";
	R[0] = 0;
	for (int i = 0; i < N; i++) {
		RR[i] = LR[i] = 0;
		RC[i] = LC[i] = 0;
		RI[i] = LI[i] = -1;
	}
	std::sort(P.begin(), P.end(), [&](const Pos& p, const Pos& q) -> bool { return p / q > 0; });
	for (int i = 0; i < N; i++) P[i].i = i, P[i].f = -1;
	Polygon H;
	for (const Pos& p : P) {
		while (H.size() > 1 && ccw(H[H.size() - 2], H.back(), p) <= 0) H.pop_back();
		H.push_back(p);
	}
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		//std::cout << "H[" << i << "]:: " << H[i] << "\n";
		P[H[i].i].f = i;
		R[i + 1] = R[i] + (H[i] - H[(i + 1) % sz]).mag();
	}
	Polygon C_;
	for (int i = 0; i < N; i++) if (P[i].f == -1) C_.push_back(P[i]);
	C_.push_back(O);
	Polygon I = graham_scan(C_);
	H.push_back(O);
	H = graham_scan(H);
	Polygon S;
	int z = 0;
	ld r = 0;
	for (int i = 0; i < N; i++) {
		if (P[i].f != -1) { r = 0; z = i; LI[i] = P[z].f; S = { P[i] }; continue; }
		LI[i] = P[z].f;
		while (S.size() > 1 && ccw(S[S.size() - 2], S.back(), P[i]) <= 0) {
			r -= (S[S.size() - 2] - S.back()).mag();
			S.pop_back();
		}
		r += (S.back() - P[i]).mag();
		LR[i] = r;
		S.push_back(P[i]);
		LC[i] = S.size() - 1;
	}
	r = 0;
	S.clear();
	for (int i = N - 1; i >= 0; i--) {
		if (P[i].f != -1) { r = 0; z = i; RI[i] = P[z].f; S = { P[i] }; continue; }
		RI[i] = P[z].f;
		while (S.size() > 1 && ccw(S[S.size() - 2], S.back(), P[i]) >= 0) {
			r -= (S[S.size() - 2] - S.back()).mag();
			S.pop_back();
		}
		r += (S.back() - P[i]).mag();
		RR[i] = r;
		S.push_back(P[i]);
		RC[i] = S.size() - 1;
	}
	r = 0;
	if (H.size() + I.size() - 2 == N) {
		int sz = H.size();
		for (int i = 0; i < sz; i++) r += (H[i] - H[(i + 1) % sz]).mag();
		sz = I.size();
		for (int i = 0; i < sz; i++) r += (I[i] - I[(i + 1) % sz]).mag();
	}
	ld l0 = P[0].mag(), l1 = P.back().mag();
	ld r0 = l0, r1 = l1;
	for (int i = 0, j = 1; i < N - 1; i++, j++) {
		r0 = l0, r1 = l1;
		r0 += LR[i];
		r0 += R[LI[i]] - R[0];
		r0 += P[i].mag();
		r1 += RR[j];
		r1 += R[sz - 1] - R[RI[j]];
		r1 += P[j].mag();
		int c = LC[i] + RC[j];
		c += LI[i] - LI[0] + 1;
		c += RI[N - 1] - RI[j] + 1;
		int cl = LC[i] + LI[i] - LI[0] + 1;
		int cr = RC[j] + RI[N - 1] - RI[j] + 1;
		if (cl > 1 && cr > 1 && c == N) r = std::max(r, r0 + r1);
	}
	return r;
}
ld no_O(Polygon& P, Polygon& H, const ld& r_) {
	std::cout << "NO::\n";
	R[0] = 0;
	for (int i = 0; i < N; i++) {
		RR[i] = LR[i] = 0;
		RC[i] = LC[i] = 0;
		RI[i] = LI[i] = -1;
	}
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		P[H[i].i].f = i;
		R[i + 1] = R[i] + (H[i] - H[(i + 1) % sz]).mag();
	}
	int s = 0;
	for (; s < N; s++) if (P[s].f != -1) break;
	Polygon S;
	int z = s;
	ld r = 0;
	for (int j = s, i = j % N; j != N + s; j++, i = j % N) {
		if (P[i].f != -1) { r = 0; z = i; LI[i] = z; S = { P[i] }; continue; }
		LI[i] = z;
		while (S.size() > 1 && ccw(S[S.size() - 2], S.back(), P[i]) <= 0) {
			r -= (S[S.size() - 2] - S.back()).mag();
			S.pop_back();
		}
		r += (S.back() - P[i]).mag();
		LR[i] = r;
		S.push_back(P[i]);
		LC[i] = S.size() - 1;
	}
	S.clear();
	z = s;
	for (int j = s, i = j % N; j != N + s; j++, i = (i - 1 + N) % N) {
		if (P[i].f != -1) { r = 0; z = i; RI[i] = z; S = { P[i] }; continue; }
		RI[i] = z;
		while (S.size() > 1 && ccw(S[S.size() - 2], S.back(), P[i]) >= 0) {
			r -= (S[S.size() - 2] - S.back()).mag();
			S.pop_back();
		}
		r += (S.back() - P[i]).mag();
		RR[i] = r;
		S.push_back(P[i]);
		RC[i] = S.size() - 1;
	}
	r = 0;
	ld r0 = 0, r1 = 0;
	auto hf = [&](const Pos& p, const Pos& q) -> bool {
		ll fc = p * q;
		ll tq = p / q;
		return !tq ? fc > 0 : tq > 0;
		};
	for (int i = 0, j = 1, i1, j1; i < N; i++) {
		while (hf(P[i], P[j])) j = (j + 1) % N;
		i1 = (j - 1 + N) % N;
		j1 = (i - 1 + N) % N;
		int s, e;
		s = RI[i]; e = LI[i1];
		r0 = s <= e ? (R[e] - R[s]) : (R[sz] - (R[s] - R[e]));
		r0 += RR[i] + LR[i1];
		r0 += P[i].mag() + P[i1].mag();
		s = RI[j]; e = LI[j1];
		r1 = s <= e ? (R[e] - R[s]) : (R[sz] - (R[s] - R[e]));
		r1 += RR[j] + LR[j1];
		r1 += P[j].mag() + P[j1].mag();
		r = std::max(r, r0 + r1);
	}
	return r;
}
void query() {
	std::cin >> N; Polygon P(N); for (Pos& p : P) std::cin >> p;
	std::sort(P.begin(), P.end(), cmpt);
	for (int i = 0; i < N; i++) P[i].i = i;
	bool f0 = 0, f1 = 1;
	Polygon C = P; C.push_back(O);
	Polygon H = graham_scan(C);
	Polygon H1 = graham_scan(C, -1);
	for (const Pos& p : H) if (p == O) { f0 = 1; break; }
	if (f0) {
		std::sort(P.begin(), P.end(), [&](const Pos& p, const Pos& q) -> bool { return p / q > 0; });
		for (int i = N - 2; i >= 0; i--) {
			if (P[N - 1] / P[i] == 0) { f0 = 0; break; }
			if (P[N - 1] / P[i]) break;
		}
	}
	f1 = H.size() == H1.size();
	bool f = 0;
	for (int i = 0, j = i + 1; i < N - 1; i++, j++)
		if (P[i] / P[j] == 0) { f = 1; break; }
	ld r = 0;
	if (f0 && f1) {
		Polygon C_ = { O };
		int sz = H.size();
		for (int i = 0; i < sz; i++) if (H[i].i != -1) P[H[i].i].f = i;
		for (int i = 0; i < N; i++) if (P[i].f == -1) C_.push_back(P[i]);
		Polygon I = graham_scan(C_);
		//std::cout << "H:: " << H.size() << " I:: " << I.size() << "\n";
		if (I.size() + H.size() - 2 == N) {
			for (int i = 0; i < sz; i++) r += (H[i] - H[(i + 1) % sz]).mag();
			int sz = I.size();
			for (int i = 0; i < sz; i++) r += (I[i] - I[(i + 1) % sz]).mag();
		}
		if (f) { std::cout << r << "\n"; return; }
	}
	std::cout << "r:: " << r << "\n";
	if (f) { std::cout << .0 << "\n"; return; }
	for (const Pos& p : H) if (p == O) { f = 1; break; }
	std::cout << (f ? yes_O(P, r) : no_O(P, H, r)) << "\n";
	return;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> Q; while (Q--) query();
	return;
}
int main() { solve(); return 0; }//boj31180

/*

1
5
1 1
2 4
1 4
0 4
-1 1

3
4
0 3
1 3
3 1
3 0
4
-4 0
5 3
0 -4
-1 0
5
4 4
5 0
3 3
3 2
-4 2
5
1 1
2 4
1 4
0 4
-1 1

2
4
-1 -1
1 1
1 -1
-1 1
5
-1 -1
1 1
1 -1
-1 1
0 2


*/

#define BOJ
#ifdef BOJ
#else
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	int Q = 130;
	for (int q = 0; q <= Q; q++) {
		std::string Din = "../../tests/aed/input/input";
		std::string Dout = "../../tests/aed/output/output";
		std::string Qin = Din + std::to_string(q) + ".txt";
		std::string Qout = Dout + std::to_string(q) + ".txt";
		freopen(Qin.c_str(), "r", stdin);
		//freopen("../../tests/aed/input/ret.txt", "w", stdout);
		Pos u;
		std::cin >> N; Polygon H(N); std::cin >> u;
		for (Pos& p : H) std::cin >> p; norm(H);
		Vint I, J;
		for (int i = 0, i1; i < N; i++) {
			i1 = (i + 1) % N;
			Pos p0 = H[i], p1 = H[i1];
			bool vis;
			if (!ccw(u, p0, p1)) {
				if (dot(p0, p1, u) < 0) std::swap(p0, p1);
				vis = 1;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					if (intersect(q0, q1, p1, u)) { vis = 0; break; }
				}
			}
			else {
				if (ccw(u, p0, p1) < 0) std::swap(p0, p1);
				Seg s1 = Seg(p0, p1);
				Vrange V;
				for (int j = 0, j1; j < N; j++) {
					if (j == i) continue;
					j1 = (j + 1) % N;
					Pos q0 = H[j], q1 = H[j1];
					if (ccw(u, q0, q1) < 0) std::swap(q0, q1);
					Seg s2 = Seg(q0, q1);
					Range r = range(s1, s2, u);
					if (r.s.den != -1) V.push_back(r);
				}
				vis = 0;
				std::sort(V.begin(), V.end());
				V.push_back(Range(Frac(1), Frac(1)));
				int sz = V.size();
				Frac hi = Frac(0);
				for (const Range& r : V) {
					if (hi < r.s) { vis = 1; break; }
					else hi = std::max(hi, r.e);
				}
			}
			if (vis) I.push_back(i + 1);
		}
		//std::cout << I.size() << "\n";
		//for (int i : I) std::cout << i << " ";
		freopen(Qout.c_str(), "r", stdin);
		std::cin >> M;
		J.resize(M); for (int& j : J) std::cin >> j;
		if (I.size() != M) { std::cout << "fuck::\n"; continue; }
		bool f = 1;
		for (int j = 0; j < M; j++) {
			if (I[j] != J[j]) { f = 0; std::cout << "fuck::\n"; break; }
		}
		if (f) std::cout << "good::\n";
		else std::cout << "fuck::\n";
	}
	return;
}
#endif
//int main() { solve(); return 0; }//boj13090
//boj 27712 10239 22635 29691 31392 16068