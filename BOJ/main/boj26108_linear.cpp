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
const ll INF = 1e17;
const int LEN = 1e5 + 1;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }
inline ll sq(const ll& x) { return x * x; }

#define TOP 0 //TOP
#define BOT 1 //BOT

#define TE 0 //TOP_E
#define TC 1 //TOP_C
#define BE 2 //BOT_E
#define BC 3 //BOT_C

#define EXC 0
#define CND 1

#define INSERT 0
#define SWAP 1

const int N_LEN = 1 << 16;
const int K_LEN = 1 << 8 | 1 << 6;

int N, K;
struct Pos {
	int x, y;
	int pi, hi, i;
	Pos(int x_ = 0, int y_ = 0, int p_ = -1, int h_ = -1, int i_ = -1) : x(x_), y(y_), pi(p_), hi(h_), i(i_) {}
	bool operator == (const Pos& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pos& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pos& p) const { return x == p.x ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const int& n) const { return { x * n, y * n }; }
	Pos operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pos& p) const { return (ll)x * p.x + (ll)y * p.y; }
	ll operator / (const Pos& p) const { return (ll)x * p.y - (ll)y * p.x; }
	Pos operator - () const { return { -x, -y }; }
	Pos operator ~ () const { return { -y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpx_rvs(const Pos& p, const Pos& q) { return p.x == q.x ? p.y > q.y : p.x < q.x; }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
Polygon monotone_chain(Polygon& C) {
	Polygon H;
	//std::sort(C.begin(), C.end());
	//C.erase(unique(C.begin(), C.end()), C.end());
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
struct Event {
	Pos v;
	int t, i, j;
	Event(Pos v_ = Pos(), int t_ = -1, int i_ = -1, int j_ = -1) : v(v_), t(t_), i(i_), j(j_) {}
	bool operator < (const Event& o) const { return v / o.v > 0; }
};
struct Jaw {
	Pos E[K_LEN];//except
	Polygon C;//candidate
	std::priority_queue<Event> EQ, CQ;
	int t, K, h, l, cnt[K_LEN];
	bool VE[N_LEN], VC[N_LEN];
	Pos ref;
	Jaw(int t_ = BOT, int k_ = -1) : t(t_), K(k_) {
		ref = (t == BOT ? Pos(0, -1) : Pos(0, 1));
		memset(VE, 0, sizeof VE);
		memset(VC, 0, sizeof VC);
		memset(cnt, 0, sizeof cnt);
		h = 0; l = 0;
	}
	void init(const Polygon& Q, const Polygon H[], const int& N, const int& K_, const int& h_, const Pos& cur = Pos(0, -1)) {
		K = K_; h = h_;
		int sz = Q.size(); assert(N > K);
		C.reserve(K_LEN);
		if (t == BOT) for (int i = 0; i <= K; i++) E[i] = Q[i];
		else for (int i = 0; i <= K; i++) {
			E[i] = Q[sz - i - 1];
			VE[E[i].pi] = 1;
			cnt[E[i].hi]++;
		}
		for (int k = 0; k <= h; k++) {
			int sz = H[k].size();
			for (int i = 0; i < sz; i++) {
				if (ccw(E[K], E[K] + cur, H[k][i]) > 0) {
					C.push_back(H[k][i]);
					VC[C[i].pi] = 1;
					l++;
					break;
				}
			}
		}
		Pos vec;
		int tq, fc;
		for (int i = 0, j; i < K; i++) {
			j = i + 1;
			vec = E[j] - E[i];
			tq = cur / vec; fc = cur * vec;
			if (tq > 0 || (!tq && fc > 0))
				EQ.push(Event(vec, (t == BOT ? BE : TE), i, j));
		}
		sz = C.size();
		for (int i = 0, j; i < sz - 1; i++) {
			j = i + 1;
			vec = C[j] - C[i];
			tq = sign(cur / vec); fc = sign(cur * vec);
			if (tq > 0 || (!tq && fc > 0))
				CQ.push(Event(vec, (t == BOT ? BC : TC), i, j));
		}
		return;
	}
	bool candidate_insert_check(const Pos& cur) {
		const Pos& e = E[K], & c = C[0];
		int tq = ccw(e, e + cur, c), fc = sign(cur * (c - e));
		return tq < 0 || (!tq && fc > 0);
	}
	bool ccw_check(const Pos& vec, const Pos& cur) {
		int tq = sign(ref / vec), fc = ref * vec;
		if (tq < 0 || (!tq && fc < 0)) return 0;
		tq = sign(cur / vec), fc = cur * vec;
		return tq > 0 || (!tq && fc > 0);
	}
	bool candidate_update(const Pos& p, const Pos& cur, int f = SWAP) {//O(K)
		Polygon tmp;
		int sz = C.size(), k = f;
		if (VE[p.i]) {
			if (f == INSERT) return 0;
			for (; k < sz; k++) tmp.push_back(C[k]);
			C = tmp;
			return 0;
		}
		VC[p.i] = 1;
		for (; k < sz; k++) {
			int tq = ccw(p, p + cur, C[k]), fc = sign(cur * (p - C[k]));
			if (tq > 0 || (!tq && fc > 0)) break;
		}
		for (int i = f; i < k; i++) tmp.push_back(C[i]);
		if (tmp.size()) {
			Pos vec = p - tmp.back();
			if (ccw_check(vec, cur)) CQ.push(Event(vec, t, k - 1, k));
		}
		tmp.push_back(p);
		if (k < sz) {
			Pos vec = p - tmp[k];
			if (ccw_check(vec, cur)) CQ.push(Event(vec, t, k, k + 1));
		}
		for (int i = k; i < sz; i++) tmp.push_back(C[i]);
		C = tmp;
		return 1;
	}
	int valid(const Event& ev, const Pos& cur, const int g = EXC) {
		Pos v;
		if (g == EXC) {
			if (!VE[E[ev.i].pi] || !VE[E[ev.j].pi]) return 2;
			v = E[ev.j] - E[ev.i];
		}
		if (g == CND) {
			if (!VC[C[ev.i].pi] || !VC[C[ev.j].pi]) return 2;
			v = C[ev.j] - C[ev.i];
		}
		int tq = sign(cur / v), fc = sign(cur * v);
		return !tq && fc > 0;
	}
	bool rotate(std::vector<Event>& EV, const Polygon H[], const Pos& cur) {
		bool f = 0;
		while (1) {
			Event ev = EQ.top();
			int val = valid(ev, cur, EXC);
			if (!val) break;
			EQ.pop();
			if (val == 2) continue;
			std::swap(E[ev.i], E[ev.j]);
			if (ev.i > 0) {
				Pos v = E[ev.i] - E[ev.i - 1];
				if (ccw_check(v, cur)) EQ.push(Event(v, t, ev.i - 1, ev.i));
			}
			if (ev.j < K) {
				Pos v = E[ev.j + 1] - E[ev.j];
				if (ccw_check(v, cur)) EQ.push(Event(v, t, ev.j, ev.j + 1));
			}
			EV.push_back(ev);
			f = 1;
		}
		while (1) {
			Event ev = CQ.top();
			int val = valid(ev, cur, CND);
			if (!val) break;
			CQ.pop();
			if (val == 2) continue;
			std::swap(C[ev.i], C[ev.j]);
			if (ev.i > 0) {
				Pos v = C[ev.i] - C[ev.i - 1];
				if (ccw_check(v, cur)) CQ.push(Event(v, t, ev.i - 1, ev.i));
			}
			if (ev.j < K) {
				Pos v = C[ev.j + 1] - C[ev.j];
				if (ccw_check(v, cur)) CQ.push(Event(v, t, ev.j, ev.j + 1));
			}
		}
		if (candidate_insert_check(cur)) {//O(NK)
			Pos bck = E[K];
			VE[bck.i] = -1;
			E[K] = C[0];
			cnt[E[K].hi]++;
			Pos v = E[K] - E[K - 1];
			if (ccw_check(v, cur)) EQ.push(Event(v, t, K - 1, K));
			int sz = H[E[K].hi].size();
			Pos p = H[E[K].hi][(E[K].i + 1) % sz];
			candidate_update(p, cur, SWAP);
			if (cnt[bck.hi] == H[bck.hi].size()) {
				candidate_update(bck, cur, INSERT);
			}
			cnt[bck.hi]--;
		}
		return f;
	}
};
struct Calipers {
	const Pos s = Pos(0, -1), e = Pos(0, 1);
	int N, K, h;
	Pos cur;
	Polygon H[K_LEN], L[K_LEN], U[K_LEN];
	Jaw bot = Jaw(BOT), top = Jaw(TOP);
	Calipers() { N = -1, K = -1, h = 0; }
	void init(Polygon& P, int n = -1, int k = -1, Pos c = Pos(0, -1)) {
		N = n; K = k; cur = c;
		std::sort(P.begin(), P.end());
		Polygon Q = P;
		Vbool F(N, 0);
		for (int k = 0; k <= K; k++) {//O(NK)
			Polygon tmp;
			for (const Pos& q : Q) if (!F[q.pi]) tmp.push_back(q);
			if (tmp.empty()) break;
			h++;
			H[k] = monotone_chain(tmp);
			int sz = H[k].size();
			for (int i = 0; i < sz; i++) H[k][i].hi = k, H[k][i].i = i;
			int b = H[k][0].i;
			for (int i = 0; i < sz; i++) if (H[k][b] < H[k][i]) b = i;
			for (int i = 0; i <= b; i++) L[k].push_back(H[k][i]);
			for (int i = b; i < sz; i++) U[k].push_back(H[k][i]);
			for (const Pos& h : H[k]) F[h.pi] = 1;
			Q.clear();
			for (const Pos& t : tmp) if (!F[t.pi]) Q.push_back(t);
		}
		Q.clear(); for (int i = 0; i < N; i++) if (F[i]) Q.push_back(P[i]);
		std::sort(Q.begin(), Q.end(), cmpx_rvs);
		bot.init(Q, L, N, K, h, s); top.init(Q, U, N, K, h, e);
	}
	bool rotate() {
		Polygon V;
		if (bot.EQ.size()) V.push_back(bot.EQ.top().v);
		if (bot.CQ.size()) V.push_back(bot.CQ.top().v);
		if (top.EQ.size()) V.push_back(-top.EQ.top().v);
		if (top.CQ.size()) V.push_back(-top.CQ.top().v);
		if (V.empty()) { cur = Pos(0, 1); return 0; }
		cur = V[0];
		int sz = V.size();
		for (int i = 1; i < sz; i++) if (cur / V[i] < 0) cur = V[i];
		return 1;
	}
	ld dist(const int& b = -1, const int& t = -1) {
		if (b < 0 || K < b || t < 0 || K < t) return INF;
		return cross(bot.E[b], bot.E[b] + cur, top.E[t]) / cur.mag();
	}
	bool jaw_rotate(ld& d) {
		int tq = e / cur;
		if (tq > 0 || (!tq && (e * cur) > 0)) return 0;
		std::vector<Event> EV;
		bot.rotate(EV, L, cur);
		top.rotate(EV, U, -cur);
		for (const Event& ev : EV) {
			if (ev.t == BOT) d = std::min(d, dist(ev.i, K - ev.i));
			else d = std::min(d, dist(K - ev.i, ev.i));
		}
		return rotate();
	}
} C;
void solve() { 
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N >> K;
	Polygon P(N); for (Pos& p : P) std::cin >> p;
	for (int i = 0; i < N; i++) P[i].i = i;
	C.init(P, N, K, Pos(0, -1));
	ld ret = INF;
	while (C.jaw_rotate(ret)) {}//O(NKlogK)
	//새로 시작하면 현재 이벤트 벡터 하나 깐다.
	//반바퀴 다 돌았으면 나간다.
	//회전에 의해 점이 새로 들어오고 나간다. 이를 처리해야 한다.
	//정렬 후 새로 들어올 가능성이 있는 점과 맨 뒤에 있는 점에 대해 현재 이벤트 벡터와의 방향성을 검사한다.
	//맨 뒤에 있는 점이 더 뒤에 있다고 판단한다면 점을 방출한 후 새로 넣는다.
	//아마도? 1대1 교환이 이루어지긴 할 것임? 새로 들어올 후보 점을 다시 한 번 검사해서 또 들어올 수 있는지를 검사해야 한다.
	//이벤트 벨트 안에 점이 정확히 K개 제외되는지 엄밀한 검사를 수행하는 로직이 필요함
	//앞대가리에 있는 각 껍질의 점들을 바로 O(1)에 조회 가능한 방법이 필요함.
	//한 번에 기록용 배열을 4개 가지고 있어야 하고, 메인 이벤트 말고 이벤트 앞대가리도 정렬이 되는 방식으로 구현을 해야할 수 있음
	//이벤트야 벡터로 저장되니까 그려러니 하는데 앞에서 들어올 후보 점을 판단하는 로직은 어떻게 단순화를 시켜야할까?
	//기울기가 같은 이벤트 벡터를 전부 힙에서 깐다.
	//이벤트 벡터에서 읽어온 순서가 변하는 두 점들을 전부 처리한다.
	//이벤트 벡터에서 읽어온 변한 순서의 점들을 대상으로 K 개 안에 있는 점들에 대해 계산한다.
	//가장 위와 아래에 있는 두 점의 이벤트 벡터 수직 벡터 사영 성분 간 길이의 절반으로 최소값을 갱신한다.
	//순서가 바뀐 두 점의 앞, 뒤를 새로 검사해서 이벤트를 새로 삽입한다. 힙에 넣으므로 자동으로 정렬된다.
	//새로 넣는 이벤트는 현재 기울기보다 무조건 정렬 순서 상 큰 기울기를 가져야 한다.
	std::cout << ret << "\n";
	return;
}
int main() { solve(); return 0; }//boj26108