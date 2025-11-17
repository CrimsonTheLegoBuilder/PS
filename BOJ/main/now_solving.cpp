#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#define right x
#define left y
typedef long long ll;
const int LEN = 2e5 + 1;
const ll MOD = 1'000'000'007;

int N, M, K, S[LEN], E[LEN];// , A[LEN];
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
};
const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
typedef std::vector<Pos> Polygon;
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { ll ret = cross(d1, d2, d3); return ret < 0 ? -1 : !!ret; }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
bool intersect(const Pos& s1, const Pos& s2, const Pos& d1, const Pos& d2, const bool f = 1) {
	bool f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0;
	bool f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0;
	if (!f) return f1 && f2;
	bool f3 = on_seg_strong(s1, s2, d1) ||
		on_seg_strong(s1, s2, d2) ||
		on_seg_strong(d1, d2, s1) ||
		on_seg_strong(d1, d2, s2);
	return (f1 && f2) || f3;
}
Pos inner_check_bi_search(const std::vector<Pos>& H, const Pos& p) {//convex
	int sz = H.size();
	//if (!sz) return Pos(-1, -1);
	//if (sz == 1) return p == H[0] ? Pos(0, 0) : Pos(-1, -1);
	//if (sz == 2) {
	//	int i1 = -1, i2 = -1;
	//	if (H[0] == p) i1 = 0;
	//	if (H[1] == p) i2 = 1;
	//	return Pos(i1, i2);
	//}
	if (cross(H[0], H[1], p) < 0 || cross(H[0], H[sz - 1], p) > 0) return Pos(-1, -1);
	//if (H[0] == p) return Pos(sz - 1, 1);
	//if (H[1] == p) return Pos(0, 2 % sz);
	//if (H[sz - 1] == p) return Pos(sz - 2, 0);
	//if (on_seg_weak(H[0], H[1], p)) return Pos(0, 1);
	//if (on_seg_weak(H[0], H[sz - 1], p)) return Pos(sz - 1, 0);
	int s = 0, e = sz - 1, m;
	while (s + 1 < e) {
		m = s + e >> 1;
		if (cross(H[0], H[m], p) >= 0) s = m;
		else e = m;
	}
	if (cross(H[s], H[e], p) < 0) return Pos(-1, -1);
	if (H[s] == p) return Pos((s - 1 + sz) % sz, e);
	if (H[e] == p) return Pos(s, (e + 1) % sz);
	if (on_seg_weak(H[s], H[e], p)) return Pos(s, e);
	return Pos(sz + 1, sz + 1);
}
Pos find_tangent_bi_search(const Polygon& H, const Pos& p) {
	int sz = H.size();
	Pos IN = Pos(sz + 1, sz + 1);
	Pos F = inner_check_bi_search(H, p);
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
			int tq = ccw(p1, cur, nxt);
			if (!f) tq *= -1;//normailze
			if (tq > 0) s = m + 1;
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
			int tq = ccw(H[0], H[k], p);
			if (!f) tq *= -1;//normailze
			if (tq > 0) s = k;
			else e = k;
		}

		//search lower hull
		int s1 = 0, e1 = s;
		while (s1 < e1) {
			m = s1 + e1 >> 1;
			const Pos& p1 = p, & cur = H[m], & nxt = H[(m + 1) % sz];
			int tq = ccw(p1, cur, nxt);
			if (!f) tq *= -1;//normailze
			if (tq > 0) s1 = m + 1;
			else e1 = m;
		}
		i1 = s1;
		if (!ccw(p, H[i1], H[(i1 + 1) % sz]) && dot(p, H[(i1 + 1) % sz], H[i1]) > 0) i1 = (i1 + 1) % sz;

		//search upper hull
		int s2 = e, e2 = sz - 1;
		while (s2 < e2) {
			m = s2 + e2 >> 1;
			const Pos& p1 = p, & cur = H[m], & nxt = H[(m + 1) % sz];
			int tq = ccw(p1, cur, nxt);
			if (!f) tq *= -1;//normailze
			if (tq < 0) s2 = m + 1;
			else e2 = m;
		}
		i2 = s2;
		if (!ccw(p, H[i2], H[(i2 + 1) % sz]) && dot(p, H[(i2 + 1) % sz], H[i2]) > 0) i2 = (i2 + 1) % sz;
	}
	if (ccw(p, H[i1], H[i2]) < 0) std::swap(i1, i2);
	return Pos(i1, i2);
}
bool inner_check(const Pos& p0, const Pos& p1, const Pos& p2, const Pos& q) { return ccw(p0, p1, q) >= 0 && ccw(p1, p2, q) >= 0 && ccw(p2, p0, q) >= 0; }
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(8);
	std::cin >> N >> M >> K;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	Polygon Q(M);
	for (int q = 0; q < M; q++) {
		std::cin >> Q[q];
		Pos tn = find_tangent_bi_search(P, Q[q]);
		S[q] = tn.right; E[q] = tn.left;
	}
	for (int k = 0, i, j; k < K; k++) {
		std::cin >> i >> j; i--; j--;
		int a = 0;
		const Pos& I = Q[i], & J = Q[j];
		//if (S[i] == -1 && S[j] == -1) { A[k] = N; continue; }
		if (S[i] == -1 && S[j] == -1) { std::cout << N << "\n"; continue; }
		if (S[j] == -1) {
			a = S[i] < E[i] ? (E[i] - S[i] + 2) : (N - (S[i] - E[i]) + 2);
			if (!ccw(I, P[S[i]], P[(S[i] + 1) % N])) a--;
			if (!ccw(I, P[E[i]], P[(E[i] - 1 + N) % N])) a--;
			//A[k] = a;
			std::cout << a << "\n";
			continue;
		}
		if (S[i] == -1) {
			a = S[j] < E[j] ? (E[j] - S[j] + 2) : (N - (S[j] - E[j]) + 2);
			if (!ccw(J, P[S[j]], P[(S[j] + 1) % N])) a--;
			if (!ccw(J, P[E[j]], P[(E[j] - 1 + N) % N])) a--;
			//A[k] = a;
			std::cout << a << "\n";
			continue;
		}
		const Pos& is = P[S[i]], & ie = P[E[i]];
		const Pos& js = P[S[j]], & je = P[E[j]];
		if (inner_check(I, is, ie, J)) {
			a = S[i] < E[i] ? (E[i] - S[i] + 2) : (N - (S[i] - E[i]) + 2);
			if (!ccw(I, P[S[i]], P[(S[i] + 1) % N])) a--;
			if (!ccw(I, P[E[i]], P[(E[i] - 1 + N) % N])) a--;
			//A[k] = a;
			std::cout << a << "\n";
			continue;
		}
		if (inner_check(J, js, je, I)) {
			a = S[j] < E[j] ? (E[j] - S[j] + 2) : (N - (S[j] - E[j]) + 2);
			if (!ccw(J, P[S[j]], P[(S[j] + 1) % N])) a--;
			if (!ccw(J, P[E[j]], P[(E[j] - 1 + N) % N])) a--;
			//A[k] = a;
			std::cout << a << "\n";
			continue;
		}
		int tq1, tq2;
		tq1 = ccw(I, J, is);
		tq2 = ccw(I, J, je);
		if (tq1 >= 0) { assert(tq2 >= 0); }
		else {
			int c = S[i] <= E[j] ? E[j] - S[i] + 1 : (N - (S[i] - E[j]) + 1);
			if (!ccw(I, P[S[i]], P[(S[i] + 1) % N])) a--;
			if (!ccw(J, P[E[j]], P[(E[j] - 1 + N) % N])) a--;
			a += c;
		}
		tq1 = ccw(J, I, js);
		tq2 = ccw(J, I, ie);
		if (tq1 >= 0) { assert(tq2 >= 0); }
		else {
			int c = S[j] <= E[i] ? E[i] - S[j] + 1 : (N - (S[j] - E[i]) + 1);
			if (!ccw(J, P[S[j]], P[(S[j] + 1) % N])) a--;
			if (!ccw(I, P[E[i]], P[(E[i] - 1 + N) % N])) a--;
			a += c;
		}
		//A[k] = a + 2;
		std::cout << a + 2 << "\n";
	}
	//for (int k = 0; k < K; k++) std::cout << A[k] << "\n";
	return;
}
int main() { solve(); return 0; }//boj34702