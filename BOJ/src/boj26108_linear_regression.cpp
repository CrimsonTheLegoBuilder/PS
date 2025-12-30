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
typedef std::vector<int> Vint;
typedef std::vector<bool> Vbool;
const ll INF = 1e17;
const int LEN = 5e4 + 1;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline bool zero(const ll& x) { return !x; }
inline ll sq(const ll& x) { return x * x; }

#define TOP 0
#define BOT 1

#define OUTL 0
#define HEAD 1
#define TAIL 2

#define INSERT 0
#define SWAP 1

#define LIMIT 8

const int N_LEN = 1 << 16;
const int K_LEN = 1 << 9;

int N, K;
struct Pos {
	ll x, y;
	int pi, hi, i;
	Pos(ll x_ = 0, ll y_ = 0, int p_ = -1, int h_ = -1, int i_ = -1) : x(x_), y(y_), pi(p_), hi(h_), i(i_) {}
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
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon P, H[K_LEN], L[K_LEN], U[K_LEN];
bool cmpx(const Pos& p, const Pos& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpx_rvs(const Pos& p, const Pos& q) { return p.x == q.x ? p.y > q.y : p.x < q.x; }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
ld vertical_dist(const Pos& p, const Pos& cur, const Pos& q) {
	if (!cur.x) return INF;
	if (ccw(p, p + cur, q) <= 0) return 0;
	ld dx = q.x - p.x;
	ld dy = dx / cur.x * cur.y;
	ld y = p.y + dy;
	return std::abs(q.y - y);
}
Polygon monotone_chain(Polygon& C) {
	Polygon H;
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
	int i, j, t;
	Event(Pos v_ = Pos(), int i_ = -1, int j_ = -1, int t_ = -1) : v(v_), i(i_), j(j_), t(t_) {}
	bool operator < (const Event& o) const { return v / o.v < 0; }
};
typedef std::priority_queue<Event> PQ;
typedef std::vector<Event> Events;
struct Jaw {
	PQ OQ, HQ, TQ, JQ, SQ;
	int t, K, h, hp, tp;
	int outl[K_LEN], head[K_LEN], tail[K_LEN], jaw[K_LEN], cnt[K_LEN];
	int ord[N_LEN][3], sts[N_LEN];
	Pos ref;
	Jaw(int t_ = BOT, int k_ = 0) : t(t_), K(k_) {
		ref = (t == BOT ? Pos(0, -1) : Pos(0, 1));
		memset(outl, -1, sizeof outl);//outlier points' idx
		memset(head, -1, sizeof head);//candidate points' idx (head)
		memset(tail, -1, sizeof tail);//candidate points' idx (tail)
		memset(jaw, -1, sizeof jaw);//rotating calipers' bot's idx
		memset(cnt, 0, sizeof cnt);//outlier points' num of each hull
		memset(ord, -1, sizeof ord);//points' order of each jaw
		memset(sts, -1, sizeof sts);//state of each point
		h = 0; hp = 0; tp = 0;
	}
	bool valid(const Pos& cur, const Pos& vec, const bool& evt = 1) {
		ll tq = cur / vec, fc = cur * vec;
		bool f1 = tq > 0 || (!tq && fc > 0);
		tq = ref / vec, fc = ref * vec;
		bool f2 = evt ? tq > 0 || (!tq && fc > 0) : 1;
		return f1 && f2;
	}
	bool nxt_check(const Pos& cur, const Pos& vec) { return cur / vec > 0 && ref / vec > 0; }
	void init(const Polygon& Q, const int& N, const int& K_, const int& h_, const Pos& cur = Pos(0, -1)) {
		K = K_; h = h_;
		int sz = Q.size(); assert(N > K);
		for (int i = 0, q; i <= K; i++) {
			q = (t == BOT ? i : sz - i - 1);
			outl[i] = Q[q].pi;
			ord[Q[q].pi][OUTL] = i;
			sts[Q[q].pi] = OUTL;
			cnt[Q[q].hi]++;
		}
		for (int k = 0; k < h; k++) jaw[k] = (t == BOT ? L[k][0].pi : U[k][0].pi);
		Polygon hd, tl;
		int tq = -1, fc = -1;
		Pos vec;
		const Pos& b = P[outl[K]];
		for (int k = 0; k < h; k++) {
			const Polygon& LH = (t == BOT ? L[k] : U[k]);
			const Polygon& UH = (t == BOT ? U[k] : L[k]);
			bool hf = 0, tf = 0;
			sz = LH.size();
			for (int i = 0; i < sz; i++) {
				if (sts[LH[i].pi] != OUTL) {
					hd.push_back(LH[i]);
					hf = 1;
					break;
				}
			}
			sz = UH.size();
			for (int i = sz - 1; i >= 0; i--) {
				if (sts[UH[i].pi] != OUTL) {
					tl.push_back(UH[i]);
					tf = 1;
					break;
				}
			}
			if (hf) hp++;
			if (tf) tp++;
			if (!hf) assert(!tf);
		}
		std::sort(hd.begin(), hd.end(), cmpx);
		if (t == TOP) std::reverse(hd.begin(), hd.end());
		sz = hd.size();
		assert(sz == hp);
		for (int i = 0; i < sz; i++) {
			head[i] = hd[i].pi;
			ord[hd[i].pi][HEAD] = i;
			assert(sts[hd[i].pi]);
			if (sts[hd[i].pi] == -1) sts[hd[i].pi] = HEAD;
			else sts[hd[i].pi] |= HEAD;
		}
		std::sort(tl.begin(), tl.end(), cmpx);
		if (t == TOP) std::reverse(tl.begin(), tl.end());
		sz = tl.size();
		assert(sz == tp);
		for (int i = 0; i < sz; i++) {
			tail[i] = tl[i].pi;
			ord[tl[i].pi][TAIL] = i;
			assert(sts[tl[i].pi]);
			if (sts[tl[i].pi] == -1) sts[tl[i].pi] = TAIL;
			else sts[tl[i].pi] |= TAIL;
		}
		for (int i = 0, j = 1; i < K; i++, j++) {
			vec = P[outl[j]] - P[outl[i]];
			if (valid(cur, vec))
				OQ.push(Event(vec, outl[i], outl[j]));
		}
		for (int i = 0, j = 1; i < hp - 1; i++, j++) {
			vec = P[head[j]] - P[head[i]];
			if (valid(cur, vec))
				HQ.push(Event(vec, head[i], head[j]));
		}
		for (int i = 0, j = 1; i < tp - 1; i++, j++) {
			vec = P[tail[j]] - P[tail[i]];
			if (valid(cur, vec))
				TQ.push(Event(vec, tail[i], tail[j]));
		}
		for (int i = 0; i < h; i++) {
			const Pos& p = P[jaw[i]];
			int sz = H[p.hi].size();
			if (sz < 2) continue;
			vec = H[p.hi][(p.i + 1) % sz] - H[p.hi][p.i];
			if (valid(cur, vec))
				JQ.push(Event(vec, H[p.hi][p.i].pi, H[p.hi][(p.i + 1) % sz].pi));
		}
		const Pos& ch = P[head[0]];
		const Pos& ct = P[tail[0]];
		vec = ch - b;
		if (valid(cur, vec))
			SQ.push(Event(vec, b.pi, ch.pi, HEAD));
		vec = ct - b;
		if (valid(cur, vec))
			SQ.push(Event(vec, b.pi, ct.pi, TAIL));
		return;
	}
	bool candidate_update(const Pos& p, const Pos& cur, PQ& CQ, int idx[], const int& f = HEAD, int o = -1, const int& pop_ = 0) {//O(K)
		int i = 0;
		int& c = (f == HEAD ? hp : tp);
		if (o != -1) {
			int x = idx[o];
			if (sts[x] == f) { sts[x] = -1; }
			else { sts[x] ^= f; }
			ord[x][f] = -1;
			for (i = o; i < c; i++) idx[i] = idx[i + 1], ord[idx[i]][f]--;
			c--;
			if (o > 0 && o < c) {
				Pos vec = P[idx[o]] - P[idx[o - 1]];
				if (valid(cur, vec)) {
					CQ.push(Event(vec, idx[o - 1], idx[o]));
				}
			}
		}
		if (pop_ || sts[p.pi] == OUTL) return 0;
		for (i = 0; i < c; i++) {
			if (valid(cur, P[idx[i]] - p, 0)) break;
		}
		for (int j = c - 1; j >= i; j--) idx[j + 1] = idx[j], ord[idx[j + 1]][f]++;
		c++;
		idx[i] = p.pi, ord[p.pi][f] = i;
		if (sts[p.pi] == -1) sts[p.pi] = f;
		else sts[p.pi] |= f;
		if (i > 0 && valid(cur, P[idx[i]] - P[idx[i - 1]])) {
			CQ.push(Event(P[idx[i]] - P[idx[i - 1]], idx[i - 1], idx[i]));
		}
		if (i < c - 1 && valid(cur, P[idx[i + 1]] - P[idx[i]])) {
			CQ.push(Event(P[idx[i + 1]] - P[idx[i]], idx[i], idx[i + 1]));
		}
		return 1;
	}
	void hull_jaw_rotate(const Pos& cur) {
		while (JQ.size()) {
			Event ev = JQ.top();
			if (cur / ev.v > 0) break;
			JQ.pop();
			if (cur / ev.v < 0) continue;
			const Pos& p = P[ev.i];
			int sz = H[p.hi].size(), i0 = (p.i + 1) % sz;
			if (sz < 2) continue;
			const Pos& p1 = H[p.hi][i0];
			if ((sts[p.pi] != -1 && (sts[p.pi] & HEAD)) && (sts[p1.pi] == -1 || sts[p1.pi] == TAIL)) {
				candidate_update(p1, cur, HQ, head, HEAD, ord[p.pi][HEAD]);
			}
			if ((sts[p.pi] != -1 && (sts[p.pi] & TAIL)) && (sts[p1.pi] == -1 || sts[p1.pi] == HEAD)) {
				candidate_update(p1, cur, TQ, tail, TAIL, ord[p.pi][TAIL]);
			}
			jaw[p.hi] = p1.pi;
			const Pos& p2 = H[p.hi][(i0 + 1) % sz];
			Pos vec = p2 - p1;
			if (valid(ref, vec)) JQ.push(Event(vec, p1.pi, p2.pi));
		}
		return;
	}
	void swap(const int& i, const int& j, int idx[], const int& f) {
		int u = ord[i][f], v = ord[j][f];
		std::swap(idx[u], idx[v]);
		ord[i][f] = v; ord[j][f] = u;
		return;
	}
	void candidate_jaw_rotate(const Pos& cur, PQ& CQ, int idx[], const int& f = HEAD) {
		Pos vec;
		while (CQ.size()) {
			Event ev = CQ.top();
			if (cur / ev.v > 0) break;
			CQ.pop();
			if (cur / ev.v < 0) continue;
			int i = ev.i, j = ev.j;
			if (sts[i] == -1 || sts[j] == -1) continue;
			if (!(sts[i] & f) || !(sts[j] & f)) continue;
			int u = ord[i][f], v = ord[j][f];
			if (u + 1 != v) continue;
			swap(i, j, idx, f);
			u = ord[i][f], v = ord[j][f];
			const int& c = (f == HEAD ? hp : tp);
			if (u < c - 1) {
				vec = P[idx[u + 1]] - P[idx[u]];
				if (valid(cur, vec)) CQ.push(Event(vec, idx[u], idx[u + 1]));
			}
			if (0 < v) {
				vec = P[idx[v]] - P[idx[v - 1]];
				if (valid(cur, vec)) CQ.push(Event(vec, idx[v - 1], idx[v]));
			}
		}
		return;
	}
	void outlier_jaw_rotate(const Pos& cur, Events& EV) {
		Pos vec;
		while (OQ.size()) {
			Event ev = OQ.top();
			if (cur / ev.v > 0) break;
			OQ.pop();
			if (cur / ev.v < 0) continue;
			int i = ev.i, j = ev.j;
			if (sts[i] != OUTL || sts[j] != OUTL) continue;
			int u = ord[i][OUTL], v = ord[j][OUTL];
			if (u + 1 != v) continue;
			EV.push_back(Event(cur, u, v, t));
			swap(i, j, outl, OUTL);
			u = ord[i][OUTL], v = ord[j][OUTL];
			if (u < K) {
				vec = P[outl[u + 1]] - P[outl[u]];
				if (valid(cur, vec)) OQ.push(Event(vec, outl[u], outl[u + 1]));
			}
			if (0 < v) {
				vec = P[outl[v]] - P[outl[v - 1]];
				if (valid(cur, vec)) OQ.push(Event(vec, outl[v - 1], outl[v]));
			}
		}
		return;
	}
	void pop(const Pos& p) {
		if (sts[p.pi] == OUTL) sts[p.pi] = -1;
		else return;
		int o = ord[p.pi][OUTL];
		assert(o == K);
		ord[p.pi][OUTL] = -1;
		outl[o] = -1;
		cnt[p.hi]--;
		return;
	}
	void push(const Pos& p) {
		ord[p.pi][OUTL] = K;
		outl[K] = p.pi;
		sts[p.pi] = OUTL;
		cnt[p.hi]++;
		return;
	}
	void outlier_candidate_swap(const Pos& cur, Events& EV) {
		Pos vec;
		vec = P[head[0]] - P[outl[K]];
		if (valid(cur, vec)) SQ.push(Event(vec, outl[K], head[0], HEAD));
		vec = P[tail[0]] - P[outl[K]];
		if (valid(cur, vec)) SQ.push(Event(vec, outl[K], tail[0], TAIL));
		while (SQ.size()) {
			Event ev = SQ.top();
			if (cur / ev.v > 0) break;
			SQ.pop();
			if (cur / ev.v < 0) continue;
			int i = ev.i, j = ev.j, typ = ev.t;
			Pos c, x = P[outl[K]];
			if (typ == HEAD) {
				if (outl[K] != i || head[0] != j) continue;
				c = P[head[0]];
			}
			else if (typ == TAIL) {
				if (outl[K] != i || tail[0] != j) continue;
				c = P[tail[0]];
			}
			EV.push_back(Event(cur, ord[x.pi][OUTL], -1, t));
			pop(x);
			int sz = H[c.hi].size();
			Pos pre = H[c.hi][(c.i - 1 + sz) % sz];
			Pos nxt = H[c.hi][(c.i + 1) % sz];
			assert(sts[c.pi] != -1);
			assert(sts[c.pi] != 0);
			if (sts[c.pi] & HEAD) {
				if (H[c.hi].size() == 1 || cnt[c.hi] + 1 == H[c.hi].size()) {
					candidate_update(c, cur, HQ, head, HEAD, ord[c.pi][HEAD], 1);//pop candi
				}
				else if (sts[nxt.pi] == OUTL) {
					if (sts[pre.pi] != OUTL && (sts[pre.pi] == -1 || !(sts[pre.pi] & HEAD)))
						candidate_update(pre, cur, HQ, head, HEAD, ord[c.pi][HEAD]);
					else assert(0);
				}
				else if (sts[nxt.pi] != OUTL) {
					if (sts[nxt.pi] == -1 || !(sts[nxt.pi] & HEAD)) candidate_update(nxt, cur, HQ, head, HEAD, ord[c.pi][HEAD]);
					else assert(0);
				}
				else {
					assert(0);
				}
			}
			if (sts[c.pi] != -1 && (sts[c.pi] & TAIL)) {
				if (H[c.hi].size() == 1 || cnt[c.hi] + 1 == H[c.hi].size()) {
					candidate_update(c, cur, TQ, tail, TAIL, ord[c.pi][TAIL], 1);//pop candi
				}
				else if (sts[pre.pi] == OUTL) {
					if (sts[nxt.pi] != OUTL && (sts[nxt.pi] == -1 || !(sts[nxt.pi] & TAIL)))
						candidate_update(nxt, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
					else assert(0);
				}
				else if (sts[pre.pi] != OUTL) {
					if (sts[pre.pi] == -1 || !(sts[pre.pi] & TAIL)) candidate_update(pre, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
					else assert(0);
				}
				else {
					assert(0);
				}
			}
			push(c);
			vec = P[outl[K]] - P[outl[K - 1]];
			if (valid(cur, vec)) {
				OQ.push(Event(vec, outl[K - 1], outl[K]));
			}
			sz = H[x.hi].size();
			pre = H[x.hi][(x.i - 1 + sz) % sz];
			nxt = H[x.hi][(x.i + 1) % sz];
			if ((sts[x.pi] == -1 || !(sts[x.pi] & TAIL)) && cnt[x.hi] + 1 == sz) {
				candidate_update(x, cur, TQ, tail, TAIL, -1);
			}
			else if (sts[pre.pi] != -1 && (sts[pre.pi] & TAIL)) {
				candidate_update(x, cur, TQ, tail, TAIL, ord[pre.pi][TAIL]);
			}
			if ((sts[x.pi] == -1 || !(sts[x.pi] & HEAD)) && cnt[x.hi] + 1 == sz) {
				candidate_update(x, cur, HQ, head, HEAD, -1);
			}
			else if (sts[nxt.pi] != -1 && (sts[nxt.pi] & HEAD)) {
				candidate_update(x, cur, HQ, head, HEAD, ord[nxt.pi][HEAD]);
			}
		}
		while (SQ.size()) SQ.pop();
		return;
	}
	void get_events(Event& o_h, Event& o_t) {
		Pos vec;
		vec = P[head[0]] - P[outl[K]];
		if (valid(ref, vec, 0)) o_h = Event(vec, outl[K], head[0], HEAD);
		else o_h = Event(Pos(0, 0));
		vec = P[tail[0]] - P[outl[K]];
		if (valid(ref, vec, 0)) o_t = Event(vec, outl[K], tail[0], TAIL);
		else o_t = Event(Pos(0, 0));
		return;
	}
	bool jaw_rotate(Events& EV, const Pos& cur) {
		candidate_jaw_rotate(cur, HQ, head, HEAD);
		candidate_jaw_rotate(cur, TQ, tail, TAIL);
		outlier_jaw_rotate(cur, EV);
		outlier_candidate_swap(cur, EV);
		hull_jaw_rotate(cur);
		return 1;
	}
};
struct Calipers {
	const Pos s = Pos(0, -1), e = Pos(0, 1);
	int N, K, h;
	Pos cur;
	Jaw bot = Jaw(BOT), top = Jaw(TOP);
	Events VB, VT;
	Calipers() { N = -1, K = -1, h = 0; }
	void init(int n = -1, int k = -1, Pos c = Pos(0, -1)) {
		N = n; K = k; cur = c;
		std::sort(P.begin(), P.end());
		for (int i = 0; i < N; i++) P[i].pi = i;
		Polygon Q = P;
		Vbool F(N, 0);
		for (int k = 0; k <= K; k++) {
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
			for (int i = b; i <= sz; i++) U[k].push_back(H[k][i % sz]);
			for (const Pos& h : H[k]) F[h.pi] = 1;
			Q.clear();
			for (const Pos& t : tmp) if (!F[t.pi]) Q.push_back(t);
		}
		std::sort(P.begin(), P.end());
		for (int k = 0; k < h; k++) {
			for (const Pos& p : H[k]) P[p.pi].hi = p.hi, P[p.pi].i = p.i;
		}
		Q.clear(); for (const Pos& p : P) if (F[p.pi]) Q.push_back(p);
		bot.init(Q, N, K, h, s); top.init(Q, N, K, h, e);
		return;
	}
	void get_events() {
		Event hev, tev;
		bot.get_events(hev, tev);
		if (hev.v != Pos(0, 0)) bot.SQ.push(hev);
		if (tev.v != Pos(0, 0)) bot.SQ.push(tev);
		top.get_events(hev, tev);
		if (hev.v != Pos(0, 0)) top.SQ.push(hev);
		if (tev.v != Pos(0, 0)) top.SQ.push(tev);
		return;
	}
	bool jaw_rotate() {
		Polygon V;
		get_events();
		if (bot.OQ.size()) V.push_back(bot.OQ.top().v);
		if (bot.JQ.size()) V.push_back(bot.JQ.top().v);
		if (bot.HQ.size()) V.push_back(bot.HQ.top().v);
		if (bot.TQ.size()) V.push_back(bot.TQ.top().v);
		if (bot.SQ.size()) V.push_back(bot.SQ.top().v);
		if (top.OQ.size()) V.push_back(-top.OQ.top().v);
		if (top.JQ.size()) V.push_back(-top.JQ.top().v);
		if (top.HQ.size()) V.push_back(-top.HQ.top().v);
		if (top.TQ.size()) V.push_back(-top.TQ.top().v);
		if (top.SQ.size()) V.push_back(-top.SQ.top().v);
		if (V.empty()) { cur = Pos(0, 1); return 0; }
		cur = V[0];
		int sz = V.size();
		for (int i = 1; i < sz; i++) if (cur / V[i] < 0) cur = V[i];
		return 1;
	}
	ld dist(const int& b = -1, const int& t = -1) {
		if (b < 0 || K < b || t < 0 || K < t) return INF;
		return vertical_dist(P[bot.outl[b]], cur, P[top.outl[t]]);
	}
	bool rotate(ld& d) {
		int tq = sign(e / cur);
		if (tq > 0 || (!tq && (e * cur) > 0)) return 0;
		std::vector<Event> EV;
		bot.jaw_rotate(EV, cur);
		top.jaw_rotate(EV, -cur);
		for (const Event& ev : EV) {
			if (ev.t == BOT) d = std::min(d, dist(ev.i, K - ev.i));
			else d = std::min(d, dist(K - ev.i, ev.i));
		}
		return jaw_rotate();
	}
} C;
ld rotating_calipers(Polygon& P) {
	std::sort(P.begin(), P.end());
	Polygon H = monotone_chain(P);
	int sz = H.size();
	if (sz < 3) return 0;
	ld ret = INF;
	for (int i = 0, j = 1; i < sz; i++) {
		while (ccw(H[i], H[(i + 1) % sz], H[j], H[(j + 1) % sz]) >= 0) {
			j = (j + 1) % sz;
		}
		Pos v = H[(i + 1) % sz] - H[i];
		ret = std::min(ret, vertical_dist(H[i], v, H[j]));
	}
	return ret * .5;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(15);
	//freopen("../../tests/G_LinearRegression/data/secret/0TRV08.in", "r", stdin);
	//freopen("../../tests/G_LinearRegression/data/secret/1CRF01.in", "r", stdin);
	//freopen("../../tests/G_LinearRegression/data/secret/2RND10.in", "r", stdin);
	//freopen("../../tests/26108/2.in", "r", stdin);
	std::cin >> N >> K;
	P.resize(N); for (Pos& p : P) std::cin >> p;
	if (N == 1) { std::cout << "0.000000000\n"; return; }
	if (N == 2 && K == 1) { std::cout << "0.000000000\n"; return; }
	if (N == 2 && K == 0) {
		if (P[0].x == P[1].x) std::cout << std::abs(P[0].y - P[1].y) * .5 << "\n";
		else std::cout << "0.000000000\n";
		return;
	}
	if (K == 0) { std::cout << rotating_calipers(P) << "\n"; return; }
	C.init(N, K, Pos(0, -1));
	ld ret = INF;
	while (C.rotate(ret)) {}
	std::cout << std::max(.0, ret * .5) << "\n";
	return;
}
int main() { solve(); return 0; }//boj26108
