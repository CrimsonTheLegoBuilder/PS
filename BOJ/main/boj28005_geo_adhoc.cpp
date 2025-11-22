#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<bool> Vbool;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 1e5 + 1;
const ld TOL = 1e-7;
inline int sign(const int& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }
inline ll sq(const ll& x) { return x * x; }
//inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
//inline bool zero(const ld& x) { return !sign(x); }
//inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
//ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }
//ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }

#define LO x
#define HI y

#define STRONG 0
#define WEAK 1

#define UP 1
#define DOWN 2

int N, M, T, Q;
struct Pos {
	int x, y;
	int i;
	Pos(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
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
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { y, x }; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
const Pos INVAL = Pos(-1, -1);
Pos ep;//point on event line
typedef std::vector<Pos> Polygon;
bool cmpt(const Pos& u, const Pos& v) {
	bool f0 = O < u;
	bool f1 = O < v;
	if (f0 != f1) return f1;
	return u / v > 0;	
}
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
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
	int i, wi, rvs;
	Seg(Pos s_ = Pos(), Pos e_ = Pos(), int i_ = -1, int wi_ = -1, int rvs_ = 0) :
		s(s_), e(e_), i(i_), wi(wi_), rvs(rvs_) {
		if (s.y > e.y) {
			std::swap(s, e);
			rvs = 1;
		}
	}
	bool operator<(const Seg& o) const {
		if (s == o.s) return e.x < o.e.x;
		if (e == o.e) return s.x < o.s.x;
		int tq1 = ccw(o.s, o.e, s);
		int tq2 = ccw(o.s, o.e, e);
		if (tq1 * tq2 >= 0) return tq1 ? tq1 > 0 : tq2 > 0;
		return ccw(s, e, o.s) < 0;
	}
	bool operator==(const Seg& o) const {return s == o.s && e == o.e; }
} seg[LEN];
typedef std::vector<Seg> Segs;
struct Trep {
	Seg l, r;
	int b, u;
	Trep(Seg l_, Seg r_, int b_ = -1, int u_ = -1) : l(l_), r(r_), b(b_), u(u_) {}
} room[LEN << 2];
class SplayTree {
	struct Node {
		Node* l;
		Node* r;
		Node* p;
		Seg key;
		Node(Seg s) : l(0), r(0), p(0) { key = s; }
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
			if (p->key == s/* eq */) break;
			if (p->key < s/* cmp */) {
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

public:
	SplayTree() : root(0) {}
	~SplayTree() { if (root) delete root; }
	void insert(Seg s) {
		if (!root) {
			root = new Node(s);
			return;
		}
		Node* p = root;
		Node** pp;
		while (1) {
			if (p->key < s/* cmp */) {
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
			int d = ccw(p->key.s, p->key.e, q);
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
	Node* find_rigtt(const Pos& q) {
		if (!root) return nullptr;
		Node* p = root;
		Node* c = nullptr;
		while (p) {
			int d = ccw(p->key.s, p->key.e, q);
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
} sp;
int D[LEN << 2], h[LEN << 2];
Vint G[LEN << 2]; int vp = 0;
int bfs(int s = 0) {
	std::queue<int> Q;
	Vbool V(vp, 0);
	Q.push(0);
	int c = 0;
	V[0] = 1;
	while (Q.size()) {
		int p = Q.back(); Q.pop();
		for (const int& v : G[p]) {
			if (!V[v]) {
				if (D[v] == DOWN) c += h[v];
				Q.push(v);
				V[v] = 1;
			}
		}
	}
	return c;
}
ll count(const Polygon& H, const int& s, const int& e, const int& n) {
	for (int i = 0; i < vp; i++) G[i].clear();
	Polygon P = { H[s] }, E;
	for (int i = e; i != s; i = (i + 1) % N) P.push_back(P[i]);
	int sz = P.size();
	for (int i = 0; i < sz; i++) P[i].i = i, E.push_back(P[i]);
	std::sort(E.begin(), E.end());
	for (int i = 0; i < sz; i++) {
		//Seg? Node?
		//해당 점의 왼쪽 변 찾기
		Segs B, U;
		for (int j = i; j < N; j++) {
			//막혀있는지 판단
			//오른쪽벽 찾기
			//아래쪽 변을 찾으면 B에 넣기
			//위쪽 변을 찾으면 U에 넣기
			//막혀있다면 break;
		}
		//아래쪽 변과 위쪽 변은 순서대로 정렬되어있는 형태이며, 아래쪽 변에는 현재 사다리꼴 방의 번호를 저장
		//위쪽 방은 현재 가장 번호가 큰 방 이후의 방 번호를 부여, 각 왼쪽 벽에 방 번호를 부여 (seg 객체를 사용)
		//가장 왼쪽 변은 -1, 가장 오른쪽 변은 INF, 나머지 공선점들은 어차피 모든 점이 정수 좌표이므로 각 방의 시작과 끝점을 정수 이벤트로 정렬 후
		//각각의 엇갈린 이벤트(방)들을 그래프로 연결
		//아래쪽 변들은 전부 sp에서 제거
		//위쪽 변들은 sp에 넣기
	}
	//그래프 완성 후 u region 세서 반환
	return bfs(0);
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(8);
	int i = 0;
	std::cin >> N; Polygon P(N); for (Pos& p : P) std::cin >> p, p.i = i, i++;
	Polygon C, F; Pos cur = O;
	for (int i = 0; i < N; i++) { C.push_back(P[(i + 1) % N] - P[i]); }
	std::sort(C.begin(), C.end(), cmpt);
	F = { cur };
	for (int i = 0; i < N; i++) {
		cur += C[i];
		while (i < N - 1 && C[i] / C[i + 1] == 0) {
			cur += C[i + 1];
			i++;
		}
		F.push_back(cur);
	}
	for (int _ = 0; _ < 2; _++) {
		int Y = -1e9, my = -1e9, diff;
		for (int i = 0; i < N; i++) P[i].i = i, Y = std::max(Y, P[i].y), my = std::max(my, F[i].y);
		Polygon H = monotone_chain(P);
		int sz = H.size(), U = 0;
		for (int i = 0; i < sz; i++) {
			if ((H[i].i + 1) % N != H[(i + 1) % sz].i) {
				int s = H[i].i, e = H[(i + 1) % sz].i;
				U += count(H, s, e, N);
			}
		}
		diff = (Y + U) - my;
		for (Pos& p : P) p = ~p;
		for (Pos& f : F) f.y += diff, f = ~f;
	}
	for (const Pos& f : F) std::cout << f << "\n";
	return;
}
int main() { solve(); return 0; }//boj28005