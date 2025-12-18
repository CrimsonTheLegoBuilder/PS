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

#define TE 0 //TOP_EXCEPT
#define TC 1 //TOP_CANDIDATE
#define TH 2 //TOP_HULL_INSERT
#define BE 3 //BOT_EXCEPT
#define BC 4 //BOT_CANDIDATE
#define BH 5 //BOT_HULL_INSERT

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
Polygon P;
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
typedef std::priority_queue<Event> PQ;
struct Jaw {
	PQ EQ, CQ, HQ;
	int t, K, h, c, bnd;
	int cnt[K_LEN], idx[K_LEN], cnd[K_LEN], ord[N_LEN];
	int hul[K_LEN];
	bool VE[N_LEN], VC[N_LEN];
	Pos ref;
	Jaw(int t_ = BOT, int k_ = -1) : t(t_), K(k_) {
		ref = (t == BOT ? Pos(0, -1) : Pos(0, 1));
		memset(VE, 0, sizeof VE);
		memset(VC, 0, sizeof VC);
		memset(cnt, 0, sizeof cnt);
		memset(idx, -1, sizeof idx);
		memset(cnd, -1, sizeof cnd);
		memset(ord, -1, sizeof ord);
		memset(hul, 0, sizeof hul);
		h = 0; c = 0; bnd = 0;
	}
	void init(const Polygon& Q, const Polygon H[], const int& N, const int& K_, const int& h_, const Pos& cur = Pos(0, -1)) {
		K = K_; h = h_;
		int sz = Q.size(); assert(N > K);
		if (t == BOT) {
			for (int i = 0; i <= K; i++) {
				idx[i] = Q[i].pi;
				ord[Q[i].pi] = i;
				cnt[Q[i].hi]++;
				hul[i] = H[i][0].i;
				//VE[E[i].pi] = 1;
			}
		}
		else {
			for (int i = 0; i <= K; i++) {
				idx[i] = Q[sz - i - 1].pi;
				ord[Q[i].pi] = i;
				cnt[Q[i].hi]++;
				hul[i] = H[i][0].i;
				//VE[E[i].pi] = 1;
			}
		}
		Polygon tmp;
		for (int k = 0; k <= h; k++) {
			int sz = H[k].size();
			bool f = 0;
			int o = -1;
			for (int i = 0; i < sz; i++) {
				if (ccw(Q[idx[K]], Q[idx[K]] + cur, H[k][i]) > 0) {
					tmp.push_back(H[k][i]);
					o = i;
					f = 1;
					//VC[C[i].pi] = 1;
					break;
				}
			}
			if (!f) c++;
			//여기서 모든 껍질의 점 하나씩을 후보로 올리고, 후보 점이 이벤트 라인 안에 있는 경우의 상한, 하한을 구한다.
			//크기가 작은 껍질이 먹혀있는 경우도 있고 완전히 벗어나있는 경우도 있다. 이런 이벤트들도 모두 고려 대상이 된다.
			//순서를 기록하는 배열에 의해 안정적으로 관리가 되고 있긴 하므로 잘 구분해서 구현하면 될 듯 함.
			//캘리퍼스 jaw는 메인, 후보군 외에 K개의 껍질에 대해 모두 돌아가고 있어야 한다.
			//후보군을 담당하는 큐와 껍질 회전을 담당하는 큐를 구분해서 사용하도록 한다.
			//변수명이 좀 길어지더라도 안 헷갈리려면 전부 뭐가 뭔지 제대로 이름을 짓기는 해야할 듯
			//메인, 앞대가리, 볼록껍질 큐 3개를 각각 관리하도록 하고
			//먹혀있는 껍질도 외부에 있는 껍질도 jaw는 회전을 하되 후보군에는 잘 구분해서 빼놓는다.
			//볼록껍질 큐는 말 그대로 볼록껍질의 jaw를 현재 이벤트 기울기보다 크도록 회전시켜준다.
			//제외되는 점군에 점을 포함시키지 않은 껍질도 따라 돌도록 하고 있다가 가장 가까운 점이 삽입이 가능한가를 알도록 한다.
			//외부에 있으면서 가장 가까이 있어서 언제든지 메인에 점을 삽입시킬 수 있는 껍질의 번호를 기억하도록 해준다.
		}
		std::sort(tmp.begin(), tmp.end(), cmpx_rvs);
		for (int i = 0; i < tmp.size(); i++) {
			cnd[i] = tmp[i].pi;
			ord[tmp[i].pi] = i;
		}
		Pos vec;
		int tq, fc;
		for (int i = 0, j; i < K; i++) {
			j = i + 1;
			vec = P[idx[j]] - P[idx[i]];
			tq = cur / vec; fc = cur * vec;
			if (tq > 0 || (!tq && fc > 0))
				EQ.push(Event(vec, (t == BOT ? BE : TE), idx[i], idx[j]));
		}
		for (int i = 0, j; i < c - 1; i++) {
			j = i + 1;
			vec = P[cnd[j]] - P[cnd[i]];
			tq = sign(cur / vec); fc = sign(cur * vec);
			if (tq > 0 || (!tq && fc > 0))
				CQ.push(Event(vec, (t == BOT ? BC : TC), cnd[i], cnd[j]));
		}
		return;
	}

	//껍질 회전 -> 후보 회전 -> 메인 삽입 이 순서로 가야 문제가 없겠지?

	bool candidate_insert_check(const Pos& cur) {
		const Pos& e = P[idx[K]], & c = P[cnd[0]];
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
		int k = f;
		if (VE[p.i]) {
			if (f == INSERT) return 0;
			for (; k < c; k++) cnd[k - 1] = cnd[k];
			c--;
			return 0;
		}
		Vint tmp;
		VC[p.i] = 1;
		for (; k < c; k++) {
			int tq = ccw(p, p + cur, P[cnd[k]]), fc = sign(cur * (p - P[cnd[k]]));
			if (tq > 0 || (!tq && fc > 0)) break;
		}
		for (int i = f; i < k; i++) tmp.push_back(cnd[i]);
		if (tmp.size()) {
			Pos vec = p - P[tmp.back()];
			if (ccw_check(vec, cur)) CQ.push(Event(vec, t, tmp.back(), p.pi));
		}
		tmp.push_back(p.pi);
		if (k < c) {
			Pos vec = P[tmp[k]] - p;
			if (ccw_check(vec, cur)) CQ.push(Event(vec, t, p.pi, tmp[k]));
		}
		for (int i = k; i < c; i++) tmp.push_back(cnd[i]);
		c = tmp.size();
		for (int i = 0; i < c; i++) cnd[i] = tmp[i], ord[tmp[i]] = i;
		return 1;
	}
	bool candidate_update_(const Pos& p, const Pos& cur, int f = -1/* f != -1 : SWAP */) {//O(K)
		const Pos& b = P[idx[K]];
		int tq = ccw(b, b + cur, p), fc = sign(dot(b, b + cur, p));
		if (tq < 0 || (!tq && fc < 0)) return 0;
		if (f == -1) {
			int i = 0;
			for (; i < c; i++) {
				tq = ccw(p, p + cur, P[cnd[i]]);
				fc = sign(dot(p, p + cur, P[cnd[i]]));
				if (tq > 0 || (!tq && fc < 0)) break;
			}
			for (int j = c; j >= i; j--) ord[cnd[j]]++, cnd[j + 1] = cnd[j];
			cnd[i] = p.pi;
			ord[p.pi] = i;
			//CQ push
		}
		return 1;
	}
	//candidate_update는 어차피 O(K) 짜리 함수니까
	//삽입과 교체 임무, 교체 시 교체되어야 할 점의 번호만 알려주면 알아서 수행하도록 만드는 게 나을 거 같음
	void hull_jaw_rotate(const Polygon H[], const Pos& cur) {
		while (1) {
			Event ev = HQ.top();
			if (cur / ev.v < 0) continue;
			if (cur / ev.v > 0) break;
			HQ.pop();
			const Pos& p = P[ev.i];
			int sz = H[p.hi].size(), i0 = (p.i + 1) % sz;
			hul[p.hi] = H[p.hi][i0].pi;
			HQ.push(Event(H[p.hi][p.i] - H[p.hi][i0], t, i0));
			candidate_update_(p, cur);
		}
	}
	int valid(const Event& ev, const Pos& cur, const int g = EXC) {
		Pos v;
		if (g == EXC) {
			if (!VE[ev.i] || !VE[ev.j]) return 2;
			v = P[ev.j] - P[ev.i];
		}
		if (g == CND) {
			if (!VC[ev.i] || !VC[ev.j]) return 2;
			v = P[ev.j] - P[ev.i];
		}
		int tq = sign(cur / v), fc = sign(cur * v);
		return !tq && fc > 0;
	}
	void swap(const int& i, const int& j, const int& t = EXC) {
		int u = ord[i], v = ord[j];
		if (t == EXC) { std::swap(idx[u], idx[v]); }
		else { std::swap(cnd[u], cnd[v]); }
		ord[i] = v; ord[j] = u;
	}
	bool rotate(std::vector<Event>& EV, const Polygon H[], const Pos& cur) {
		bool f = 0;
		while (1) {
			Event ev = EQ.top();
			int val = valid(ev, cur, EXC);
			if (!val) break;
			EQ.pop();
			if (val == 2) continue;
			swap(ev.i, ev.j);
			if (ord[ev.i] > 0) {
				Pos v = P[ev.i] - P[idx[ord[ev.i] - 1]];
				if (ccw_check(v, cur)) EQ.push(Event(v, t, idx[ord[ev.i] - 1], ev.i));
			}
			if (ord[ev.j] < K) {
				Pos v = P[idx[ord[ev.j] + 1]] - P[ev.j];
				if (ccw_check(v, cur)) EQ.push(Event(v, t, ev.j, idx[ord[ev.j] + 1]));
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
			swap(ev.i, ev.j);
			if (ord[ev.i] > 0) {
				Pos v = P[ev.i] - P[cnd[ord[ev.i] - 1]];
				if (ccw_check(v, cur)) CQ.push(Event(v, t, cnd[ord[ev.i] - 1], ev.i));
			}
			if (ord[ev.j] < c - 1) {
				Pos v = P[cnd[ord[ev.j] + 1]] - P[ev.j];
				if (ccw_check(v, cur)) CQ.push(Event(v, t, ev.j, cnd[ord[ev.j] + 1]));
			}
		}
		if (candidate_insert_check(cur)) {//O(NK)
			Pos bck = P[idx[K]];
			VE[idx[K]] = -1;
			idx[K] = cnd[0];
			cnt[P[idx[K]].hi]++;
			Pos v = P[idx[K]] - P[idx[K - 1]];
			if (ccw_check(v, cur)) EQ.push(Event(v, t, idx[K - 1], idx[K]));
			int sz = H[P[idx[K]].hi].size();
			Pos p = H[P[idx[K]].hi][(P[idx[K]].i + 1) % sz];
			candidate_update(p, cur, SWAP);
			if (cnt[bck.hi] == H[bck.hi].size()) candidate_update(bck, cur, INSERT);
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
	void init(int n = -1, int k = -1, Pos c = Pos(0, -1)) {
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
		std::sort(P.begin(), P.end(), cmpx_rvs);
		for (int k = 0; k < h; k++) {
			for (const Pos& p : H[k]) P[p.pi].hi = p.hi, P[p.pi].i = p.i;
		}
		return;
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
		return cross(P[b], P[b] + cur, P[t]) / cur.mag();
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
ld rotating_calipers(Polygon& P) {
	Polygon H = monotone_chain(P);
	int sz = H.size();
	ld ret = INF;
	for (int i = 0, j = 1; i < sz; i++) {
		while (ccw(H[i], H[(i + 1) % sz], H[j], H[(j + 1) % sz]) >= 0) {
			j = (j + 1) % sz;
		}
		Pos v = H[(i + 1) % sz] - H[i];
		ld tq = cross(H[i], H[(i + 1) % sz], H[j]);
		ret = std::min(ret, tq / v.mag());
	}
	return ret;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	std::cin >> N >> K;
	P.resize(N); for (Pos& p : P) std::cin >> p;
	if (N < 3) { std::cout << "0.000000000\n"; return; }
	if (N <= 3 && K == 1) { std::cout << "0.000000000\n"; return; }
	if (K == 0) { std::cout << rotating_calipers(P) << "\n"; return; }
	std::sort(P.begin(), P.end(), cmpx_rvs);
	for (int i = 0; i < N; i++) P[i].i = i;
	C.init(N, K, Pos(0, -1));
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


	//왜 볼록껍질들에 대해서도 관리가 필요한가 했는데
	//처음에는 이벤트 구역에 걸쳐있지 않다가 회전 중에 걸치는 경우가 생긴다.
	//구조체 재정의는 할 필요는 없을 거 같기는 한데... 이미 타입을 나타내는 멤버변수가 있다.
	//문제는 어떻게 구현을 하는가 이다. jaw 안에 이중jaw를 구현해야한다.
	//후보군의 정렬 순서를 기억하고 있게 하되
	//완전히 먹힌 경우와
	//아예 벗어나있는 경우를 각각 구분할 수 있도록 하여야 함.
	//완전히 먹힌 경우는 현재 이벤트 벨트에 들어가있는 점의 수를 세기 때문에 원천 배제가 가능하다고 치고
	//외부에 모든 껍질이 드러나있는 경우는 정렬 후 이벤트를 계속 돌려야 함.
	//완전히 먹힌 경우도 마찬가지로 이벤트는 계속 만들어서 돌려야 함.

	std::cout << ret << "\n";
	return;
}
int main() { solve(); return 0; }//boj26108