#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
typedef long long ll;
//typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<bool> Vbool;
const ll INF = 1e9;
const int LEN = 1e5 + 1;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }

#define UP 1
#define DOWN 2

int N;
struct Pos {
	ll x, y;//fucking overflow
	int i;
	Pos(ll x_ = 0, ll y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos operator ! () const { return { y, x }; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
bool cmpy(const Pos& p, const Pos& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpt(const Pos& u, const Pos& v) {
	bool f0 = O < u;
	bool f1 = O < v;
	if (f0 != f1) return f1;
	return u / v > 0;
}
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && dot(d1, d3, d2) > 0; }
Polygon monotone_chain(Polygon C) {
	Polygon H;
	std::sort(C.begin(), C.end());
	C.erase(unique(C.begin(), C.end()), C.end());
	if (C.size() <= 2) { for (const Pos& pos : C) H.push_back(pos); }
	else {
		for (int i = 0; i < C.size(); i++) {
			while (H.size() > 1 && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) < 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
		int s = H.size() + 1;
		for (int i = C.size() - 1; i >= 0; i--) {
			while (H.size() > s && ccw(H[H.size() - 2], H[H.size() - 1], C[i]) < 0)
				H.pop_back();
			H.push_back(C[i]);
		}
		H.pop_back();
	}
	return H;
}
struct Seg {
	Pos s, e;
	int i, ri, rvs;
	Seg(Pos s_ = Pos(), Pos e_ = Pos(), int i_ = -1, int ri_ = -1, int rvs_ = 0) :
		s(s_), e(e_), i(i_), ri(ri_), rvs(rvs_) {
		if (s.y > e.y) { std::swap(s, e); rvs = 1; }
	}
	bool operator < (const Seg& o) const {
		if (s == o.s) return ccw(o.s, o.e, e) > 0;
		if (e == o.e) return ccw(o.s, o.e, s) > 0;
		int tq1 = ccw(o.s, o.e, s);
		int tq2 = ccw(o.s, o.e, e);
		if (tq1 * tq2 >= 0) return tq1 ? tq1 > 0 : tq2 > 0;
		return ccw(s, e, o.s) < 0;
	}
	bool operator == (const Seg& o) const { return s == o.s && e == o.e; }
} seg[LEN];
typedef std::vector<Seg> Segs;
bool on_seg_strong(const Seg& se, const Pos& q) { return !ccw(se.s, se.e, q) && dot(se.s, q, se.e) >= 0; }
struct Trap {
#ifdef STRICT
	int l, r;
	ll b, u;
	Trap(int l_ = -1, int r_ = -1, ll b_ = -1, ll u_ = -1) : l(l_), r(r_), b(b_), u(u_) {}
#else
	ll b, u;
	Trap(ll b_ = -1, ll u_ = -1) : b(b_), u(u_) {}
#endif
	ll h() const { return u - b; }
} room[LEN << 2];
class SplayTree {
	struct Node {
		Node* l;
		Node* r;
		Node* p;
		int val;
		Node(int i) : l(0), r(0), p(0) { val = i; }
		~Node() { if (l) delete l; if (r) delete r; }
	} *root;
	void rotate(Node* x) {
		Node* p = x->p;
		if (!p) return;
		Node* b = 0;
		if (x == p->l) {
			p->l = b = x->r;
			x->r = p;
		}
		else {
			p->r = b = x->l;
			x->l = p;
		}
		x->p = p->p;
		p->p = x;
		if (b) b->p = p;
		(x->p ? p == x->p->l ? x->p->l : x->p->r : root) = x;
	}
	void splay(Node* x) {
		while (x->p) {
			Node* p = x->p;
			Node* g = p->p;
			if (g) {
				if ((x == p->l) == (p == g->l)) rotate(p);
				else rotate(x);
			}
			rotate(x);
		}
	}
	Node* find(Seg s) {
		if (!root) return 0;
		Node* p = root;
		while (1) {
			if (seg[p->val] == s) break;
			if (seg[p->val] < s) {
				if (!p->r) {
					break;
				}
				p = p->r;
			}
			else {
				if (!p->l) {
					break;
				}
				p = p->l;
			}
		}
		splay(p);
		return p;
	}

private:
	void _delete_recursive(Node* node) {
		if (!node) return;
		_delete_recursive(node->l);
		_delete_recursive(node->r);
		delete node;
	}

public:
	SplayTree() : root(0) {}
	~SplayTree() { if (root) delete root; }
	void clear() { if (root) { _delete_recursive(root); root = nullptr; } }
	bool empty() const { return root == nullptr; }
	void insert(int s) {
		if (!root) {
			root = new Node(s);
			return;
		}
		Node* p = root;
		Node** pp;
		while (1) {
			if (seg[p->val] < seg[s]) {
				if (!p->r) {
					pp = &p->r;
					break;
				}
				p = p->r;
			}
			else {
				if (!p->l) {
					pp = &p->l;
					break;
				}
				p = p->l;
			}
		}
		Node* x = new Node(s);
		*pp = x;
		x->p = p;
		splay(x);
	}
	void pop(Node* ptr) {
		if (!ptr) return;
		splay(ptr);
		Node* p = root;
		if (p->l && p->r) {
			root = p->l; root->p = 0;
			Node* l = root;
			while (l->r) l = l->r;
			l->r = p->r;
			p->r->p = l;
		}
		else if (p->l) root = p->l, root->p = 0;
		else if (p->r) root = p->r, root->p = 0;
		else root = 0;
		p->l = p->r = 0;
		delete p;
	}
	int erase(int i) {
		if (!find(seg[i])) return 0;
		Node* p = root;
		if (p->l && p->r) {
			root = p->l; root->p = 0;
			Node* l = root;
			while (l->r) l = l->r;
			l->r = p->r;
			p->r->p = l;
		}
		else if (p->l) root = p->l, root->p = 0;
		else if (p->r) root = p->r, root->p = 0;
		else root = 0;
		p->l = p->r = 0;
		delete p;
		return 1;
	}
	Node* get_prev(Node* x) {
		if (!x || !x->l) return nullptr;
		Node* p = x->l;
		while (p->r) p = p->r;
		return p;
	}
	Node* get_next(Node* x) {
		if (!x || !x->r) return nullptr;
		Node* p = x->r;
		while (p->l) p = p->l;
		return p;
	}
	Node* find_left(const Pos& q) {
		if (!root) return nullptr;
		Node* p = root;
		Node* c = nullptr;
		while (p) {
			int d = ccw(seg[p->val].s, seg[p->val].e, q);
			if (d < 0) {
				c = p;
				p = p->r;
			}
			else {
				p = p->l;
			}
		}
		if (c) splay(c);
		return c;
	}
	bool find_left(const Pos& q, int& s) {
		Node* c = find_left(q);
		if (c) { s = c->val; return 1; }
		return 0;
	}
	Node* find_right(const Pos& q) {
		if (!root) return nullptr;
		Node* p = root;
		Node* c = nullptr;
		while (p) {
			int d = ccw(seg[p->val].s, seg[p->val].e, q);
			if (d > 0) {
				c = p;
				p = p->l;
			}
			else {
				p = p->r;
			}
		}
		if (c) splay(c);
		return c;
	}
	bool find_right(const Pos& q, int& s) {
		Node* c = find_right(q);
		if (c) { s = c->val; return 1; }
		return 0;
	}
} sp;
struct Info { int i, d; };
ll h[LEN << 2];
std::vector<Info> G[LEN << 2]; int vp = 0;
ll bfs(int s = 0) {
	std::queue<int> Q;
	Vbool V(vp, 0);
	Q.push(s);
	ll c = 0;
	V[s] = 1;
	while (Q.size()) {
		int p = Q.front(); Q.pop();
		for (const Info& w : G[p]) {
			if (!V[w.i]) {
				if (w.d == DOWN) c += h[w.i];
				Q.push(w.i);
				V[w.i] = 1;
			}
		}
	}
	return c;
}
bool block(const Polygon& P, const int& i, const int& sz) {
	int i0 = (i - 2 + sz) % sz, i1 = (i - 1 + sz) % sz, i2 = i, i3 = (i + 1) % sz, i4 = (i + 2) % sz;
	if (P[i1].y < P[i2].y && P[i2].y < P[i3].y) return 1;
	if (P[i1].y == P[i2].y && P[i1].x < P[i2].x && P[i0].y < P[i1].y && P[i2].y < P[i3].y) return 1;
	if (P[i2].y == P[i3].y && P[i3].x < P[i2].x && P[i1].y < P[i2].y && P[i3].y < P[i4].y) return 1;
	return 0;
}
struct Proj { ll s, e; int i; };
ll count(const Polygon& H, const int& s, const int& e, const int& n) {
	sp.clear();
	for (int i = 0; i < vp; i++) G[i].clear();
	vp = 1;
	Polygon P = { H[s] }, E;
	for (int i = e; i != s; i = (i - 1 + N) % N) P.push_back(H[i]);
	int sz = P.size();
	for (int i = 0; i < sz; i++) P[i].i = i, E.push_back(P[i]);
	for (int i = 0; i < sz; i++) seg[i] = Seg(P[i], P[(i + 1) % sz], i);
	std::sort(E.begin(), E.end(), cmpy);
	for (int i = 0, j, i0, i1, i2, i3, i4, y; i < sz; i++) {
		int l = -1, r = -1;
		Pos cur = E[i];
		y = cur.y;
		i2 = cur.i;
		i0 = (i2 - 2 + sz) % sz; i1 = (i2 - 1 + sz) % sz; i3 = (i2 + 1) % sz; i4 = (i2 + 2) % sz;
		if (sp.empty() ||
			(P[i1].y > y && P[i3].y > y && ccw(P[i1], P[i2], P[i3]) > 0) ||
			(P[i2].y == P[i3].y && P[i1].y > y && P[i4].y > y && ccw(P[i1], P[i2], P[i3]) > 0)
			) {
			l = i1; seg[i1] = Seg(P[i1], P[i2], i1, vp);
			if (P[i2].y == P[i3].y) r = i3, seg[i3] = Seg(P[i3], P[i4], i3, vp), i++;
			else r = i2, seg[i2] = Seg(P[i2], P[i3], i2, vp);
#ifdef STRICT
			room[vp] = Trap(l, r, y);
#else
			room[vp] = Trap(y);
#endif
			if (P[i2].y == P[i3].y && P[i2].i == 0) { G[0].push_back({ seg[l].ri, UP }); h[seg[l].ri] = 0; }
			sp.insert(l); sp.insert(r);
			vp++;
			continue;
		}
		if ((P[i1].y < y && P[i3].y < y && ccw(P[i1], P[i2], P[i3]) > 0) ||
			(P[i1].y == P[i2].y && P[i0].y < y && P[i3].y < y && ccw(P[i1], P[i2], P[i3]) > 0)
			) {
			l = i2;
			if (P[i1].y == P[i2].y) r = i0, i++;
			else r = i1;
			room[seg[l].ri].u = y;
			h[seg[l].ri] = room[seg[l].ri].h();
			if (!seg[l].i || !seg[r].i) { G[0].push_back({ seg[l].ri, UP }); h[seg[l].ri] = 0; }
			if (P[i1].y == P[i2].y && P[i1].i == 0) { G[0].push_back({ seg[l].ri, DOWN }); }
			sp.erase(l); sp.erase(r);
			continue;
		}
		bool fl = sp.find_left(E[i], l);
		bool fr = sp.find_right(E[i], r);
		if (fl && !seg[l].rvs) fl = 0;
		if (fr && seg[r].rvs) fr = 0;
		Vint B = {}, U = {};
		int tog = 0;
		if (fl && seg[l].s.y < y && y < seg[l].e.y) { B.push_back(l), U.push_back(l); }
		for (j = i; j < sz; j++) {
			if (E[j].y != y) break;
			if (fr && seg[r].s.y < y && y < seg[r].e.y && ccw(seg[r].s, seg[r].e, E[j]) < 0) break;
			cur = E[j]; i = j;
			fr = sp.find_right(cur, r);
			if (fr && seg[r].rvs) fr = 0;
			bool blk = block(P, cur.i, sz);
			int x1 = cur.i, x0 = (x1 - 1 + sz) % sz, x2 = (x1 + 1) % sz;
			int s0 = sign(P[x0].y - P[x1].y), s1 = sign(P[x2].y - P[x1].y);
			if (s0 > 0 && s1 > 0) U.push_back(x1), U.push_back(x0);
			else if (s0 < 0 && s1 < 0) B.push_back(x0), B.push_back(x1);
			else {
				if (s0 > 0) U.push_back(x0);
				else if (s0 < 0) B.push_back(x0);
				if (s1 > 0) U.push_back(x1);
				else if (s1 < 0) B.push_back(x1);
			}
			if (blk) { tog = 1; break; }
		}
		if (!tog) {
			fr = sp.find_right(cur, r);// assert(fr && seg[r].s.y < y && y < seg[r].e.y);
			B.push_back(r); U.push_back(r);
		}
		int szu = U.size(), szb = B.size();
		//assert(szu % 2 == 0);
		//assert(szb % 2 == 0);
		std::vector<Proj> PU, PB;
		for (int j = 0; j < szb; j += 2) {
			PB.push_back({ seg[B[j]].e.x, seg[B[j + 1]].e.x, seg[B[j]].ri });
			room[seg[B[j]].ri].u = y;
			h[seg[B[j]].ri] = room[seg[B[j]].ri].h();
			if (!seg[B[j]].i || !seg[B[j + 1]].i) { G[0].push_back({ seg[B[j]].ri, UP }); h[seg[B[j]].ri] = 0; }
			//assert(seg[B[j]].rvs && !seg[B[j + 1]].rvs);
		}
		for (int j = 0; j < szu; j += 2) {
			PU.push_back({ seg[U[j]].s.x, seg[U[j + 1]].s.x, vp });
			seg[U[j]].ri = vp;
			seg[U[j + 1]].ri = vp;
#ifdef STRICT
			room[vp] = Trap(U[j], U[j + 1], y);
#else
			room[vp] = Trap(y);
#endif
			vp++;
			//assert(seg[U[j]].rvs && !seg[U[j + 1]].rvs);
		}
		//assert(U.size() && B.size());
		if (U[0] == B[0]) PU[0].s = PB[0].s = -INF;
		if (U.back() == B.back()) PU.back().e = PB.back().e = INF;
		int idx_u = 0, idx_b = 0;
		while (idx_u < PU.size() && idx_b < PB.size()) {
			const Proj& u = PU[idx_u];
			const Proj& b = PB[idx_b];
			if (std::max(u.s, b.s) < std::min(u.e, b.e)) {
				G[b.i].push_back({ u.i, UP });
				G[u.i].push_back({ b.i, DOWN });
			}
			if (u.e < b.e) idx_u++;
			else idx_b++;
		}
		for (const int& b : B) sp.erase(b);
		for (const int& u : U) sp.insert(u);
	}
	return bfs();
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(8);
	std::cin >> N; Polygon P(N), H, F; for (Pos& p : P) std::cin >> p;
	Vbool V(N, 1);
	Pos cur = O;
	for (int i = 0; i < N; i++) {
		const Pos& p0 = P[(i - 1 + N) % N], & p1 = P[i], & p2 = P[(i + 1) % N];
		if (on_seg_weak(p0, p2, p1)) V[i] = 0;
	}
	for (int i = 0; i < N; i++) if (V[i]) H.push_back(P[i]);
	P = H; N = P.size(); H.clear();
	for (int i = 0; i < N; i++) { H.push_back(P[(i + 1) % N] - P[i]); }
	std::sort(H.begin(), H.end(), cmpt);
	F = { cur };
	for (int i = 0; i < N; i++) {
		cur += H[i];
		while (i < N - 1 && H[i] / H[i + 1] == 0) {
			cur += H[i + 1];
			i++;
		}
		F.push_back(cur);
	}
	F.pop_back();
	auto min_it = std::min_element(F.begin(), F.end());
	std::rotate(F.begin(), min_it, F.end());
	for (int _ = 0, sz; _ < 2; _++) {
		ll Y = -1e9, my = -1e9, diff;
		for (int i = 0; i < N; i++) P[i].i = i, Y = std::max(Y, P[i].y);
		sz = F.size();
		for (int i = 0; i < F.size(); i++) my = std::max(my, F[i].y);
		H = monotone_chain(P);
		ll U = 0;
		sz = H.size();
		for (int i = 0; i < sz; i++) {
			if ((H[i].i + 1) % N != H[(i + 1) % sz].i) {
				U += count(P, H[i].i, H[(i + 1) % sz].i, N);
			}
		}
		diff = (Y + U) - my;
		for (Pos& p : P) p = !p;
		for (Pos& f : F) f.y += diff, f = !f;
		std::reverse(P.begin(), P.end());
	}
	std::cout << F.size() << "\n";
	for (const Pos& f : F) std::cout << f << "\n";
	return;
}
int main() { solve(); return 0; }//boj28005
