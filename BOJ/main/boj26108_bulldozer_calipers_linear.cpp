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

#define IGN (-1)
#define OUTL 0
#define HEAD 1
#define TAIL 2
#define BOTH 3

#define INSERT 0
#define SWAP 1

#define LIMIT 8

const int N_LEN = 1 << 16;
const int K_LEN = 1 << 9;

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
	//std::cout << "DEBUG::\nDEBUG::\nDEBUG::\n";
	//std::cout << "	vertical dist:: p  :: " << p << "\n";
	//std::cout << "	vertical dist:: q  :: " << q << "\n";
	//std::cout << "	vertical dist:: cur:: " << cur << "\n";
	//std::cout << "DEBUG::\nDEBUG::\nDEBUG::\n";
	if (!cur.x) return INF;
	if (ccw(p, p + cur, q) <= 0) return 0;
	ld dx = q.x - p.x;
	ld dy = dx / cur.x * cur.y;
	ld y = p.y + dy;
	std::cout << "DEBUG::\nDEBUG::\nDEBUG::\n";
	std::cout << "vertical dist:: dx  :: " << dx << "\n";
	std::cout << "vertical dist:: dy  :: " << dy << "\n";
	std::cout << "vertical dist:: cur :: " << cur << "\n";
	std::cout << "vertical dist:: y   :: " << y << "\n";
	std::cout << "vertical dist:: diff:: " << q.y - y << "\n";
	std::cout << "DEBUG::\nDEBUG::\nDEBUG::\n";
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
void print_event(const Event& e) {
	std::cout << "	[v:(" << e.v.x << "," << e.v.y << ") " << "i:" << e.i << " P[i]: (" << P[e.i].x << "," << P[e.i].y << ") j:" << e.j << " P[j]: (" << P[e.j].x << "," << P[e.j].y << ") t:" << e.t << "]";
}
void debug_pq(PQ pq, const char* name = "PQ") {
	std::cout << "--- DEBUG:: [" << name << "] (Size: " << pq.size() << ") ---\n";
	int idx = 0;
	while (!pq.empty()) {
		Event e = pq.top();
		pq.pop();
		std::cout << "#" << idx++ << ": ";
		print_event(e);
		std::cout << "\n";
	}
	std::cout << "--- DEBUG:: ------------------------\n";
}
typedef std::vector<Event> Events;
struct Jaw {
	PQ OQ, HQ, TQ, JQ, SQ;
	int t, K, h, hp, tp, bnd;
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
		h = 0; hp = 0; tp = 0; bnd = 0;
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
		K = K_; h = h_;//bnd = 0;
		int sz = Q.size(); assert(N > K);
		//std::cout << "jaw init start::\n";
		for (int i = 0, q; i <= K; i++) {
			q = (t == BOT ? i : sz - i - 1);
			outl[i] = Q[q].pi;
			ord[Q[q].pi][OUTL] = i;
			sts[Q[q].pi] = OUTL;
			cnt[Q[q].hi]++;
			//bnd = std::max(bnd, Q[q].hi);
		}
		//std::cout << "outl done::\n";
		//std::cout << "h:: " << h << "\n";
		for (int k = 0; k < h; k++) jaw[k] = (t == BOT ? L[k][0].pi : U[k][0].pi);
		//std::cout << "jaw numbering done::\n";
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
				//Pos vec = LH[i] - b;
				//if (valid(cur, vec)) {
				//	hd.push_back(LH[i]);
				//	hf = 1;
				//	break;
				//}
				if (sts[LH[i].pi] != OUTL) {
					hd.push_back(LH[i]);
					hf = 1;
					break;
				}
			}
			sz = UH.size();
			for (int i = sz - 1; i >= 0; i--) {
				//Pos vec = UH[i] - b;
				//if (valid(cur, vec)) {
				//	tl.push_back(UH[i]);
				//	tf = 1;
				//	break;
				//}
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
		//std::cout << "head, tail get done::\n";
		std::sort(hd.begin(), hd.end(), cmpx_rvs);
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
		//std::cout << "head numbering done::\n";
		std::sort(tl.begin(), tl.end(), cmpx_rvs);
		if (t == TOP) std::reverse(tl.begin(), tl.end());
		//std::cout << "DEBUG:: tl:: \n";
		//for (const Pos& tt : tl) std::cout << "	" << tt << "\n";
		//std::cout << "DEBUG:: tl:: \n";
		sz = tl.size();
		assert(sz == tp);
		for (int i = 0; i < sz; i++) {
			tail[i] = tl[i].pi;
			ord[tl[i].pi][TAIL] = i;
			assert(sts[tl[i].pi]);
			if (sts[tl[i].pi] == -1) sts[tl[i].pi] = TAIL;
			else sts[tl[i].pi] |= TAIL;
		}
		//std::cout << "tail numbering done::\n";
		for (int i = 0, j = 1; i < K; i++, j++) {
			vec = P[outl[j]] - P[outl[i]];
			if (valid(cur, vec))
				OQ.push(Event(vec, outl[i], outl[j]));
		}
		//std::cout << "OQ init done::\n";
		for (int i = 0, j = 1; i < hp - 1; i++, j++) {
			vec = P[head[j]] - P[head[i]];
			if (valid(cur, vec))
				HQ.push(Event(vec, head[i], head[j]));
		}
		//std::cout << "HQ init done::\n";
		for (int i = 0, j = 1; i < tp - 1; i++, j++) {
			vec = P[tail[j]] - P[tail[i]];
			if (valid(cur, vec))
				TQ.push(Event(vec, tail[i], tail[j]));
		}
		//std::cout << "TQ init done::\n";
		for (int i = 0; i < h; i++) {
			const Pos& p = P[jaw[i]];
			//std::cout << "p.hi:: " << p.hi << "\n";
			int sz = H[p.hi].size();
			vec = H[p.hi][(p.i + 1) % sz] - H[p.hi][p.i];
			if (valid(cur, vec))
				JQ.push(Event(vec, H[p.hi][p.i].pi, H[p.hi][(p.i + 1) % sz].pi));
		}
		//std::cout << "JQ init done::\n";
		//const Pos& b = P[outl[K]];
		const Pos& ch = P[head[0]];
		const Pos& ct = P[tail[0]];
		vec = ch - b;
		if (valid(cur, vec))
			SQ.push(Event(vec, b.pi, ch.pi, HEAD)); 
		vec = ct - b;
		if (valid(cur, vec))
			SQ.push(Event(vec, b.pi, ct.pi, TAIL));
		//std::cout << (t == BOT ? "BOT::\n" : "TOP::\n");
		//debug_pq(OQ, "OQ");
		//debug_pq(HQ, "HQ");
		//debug_pq(TQ, "TQ");
		//debug_pq(JQ, "JQ");
		//debug_pq(SQ, "SQ");
		return;
	}
	bool candidate_update(const Pos& p, const Pos& cur, PQ& CQ, int idx[], const int& f = HEAD, int o = -1, const int& pop_ = 0) {//O(K)
		std::cout << "candidate update::\n";
		//std::cout << "	outl[K]:: " << outl[K] << "\n";
		//const Pos& b = P[outl[K]];
		int i = 0;
		int& c = (f == HEAD ? hp : tp);
		if (o != -1) {
			int x = idx[o];
			if (sts[x] == f) { sts[x] = -1; }
			else { sts[x] ^= f; }
			ord[x][f] = -1;
			for (i = o; i < c; i++) idx[i] = idx[i + 1], ord[idx[i]][f]--;
			c--;
		}
		std::cout << "candi update:: c:: " << c << "\n";
		std::cout << "candi update:: o:: " << o << "\n";
		for (int j = 0; j < c; j++) {
			std::cout << (f == HEAD ? "head[" : "tail[") << j << "]:: " << idx[j] << " P[idx[" << j << "]]:: " << P[idx[j]] << "\n";
		}
		//int tq = ccw(b, b + cur, p), fc = sign(dot(b, b + cur, p));
		//if (pp || !valid(cur, p - b) || sts[p.pi] == OUTL) return 0;
		if (pop_ || sts[p.pi] == OUTL) return 0;
		for (i = 0; i < c; i++) {
			std::cout << "	DEBUG:: idx[i]:: " << idx[i] << "\n";
			std::cout << "	DEBUG:: cur   :: " << cur << "\n";
			std::cout << "	DEBUG:: vec   :: " << P[idx[i]] - p << "\n";
			std::cout << "	DEBUG:: valid :: " << valid(cur, P[idx[i]] - p, 0) << "\n";
			if (valid(cur, P[idx[i]] - p, 0)) break;
		}
		std::cout << "candi update:: i:: " << i << "\n";
		for (int j = c - 1; j >= i; j--) idx[j + 1] = idx[j], ord[idx[j + 1]][f]++;
		c++;
		//std::cout << "FUCK::\n";
		idx[i] = p.pi, ord[p.pi][f] = i;
		//std::cout << "DEBUG:: p  :: " << p << "\n";
		//std::cout << "DEBUG:: idx[i]  :: " << idx[i] << "\n";
		if (sts[p.pi] == -1) sts[p.pi] = f;
		else sts[p.pi] |= f;
		if (i > 0 && valid(cur, P[idx[i]] - P[idx[i - 1]])) {
			//std::cout << "	DEBUG:: idx[i-1]:: " << idx[i - 1] << "\n";
			//std::cout << "	DEBUG:: idx[i]  :: " << idx[i] << "\n";
			CQ.push(Event(P[idx[i]] - P[idx[i - 1]], idx[i - 1], idx[i]));
		}
		if (i < c - 1 && valid(cur, P[idx[i + 1]] - P[idx[i]])) {
			//std::cout << "	DEBUG:: idx[i]  :: " << idx[i] << "\n";
			//std::cout << "	DEBUG:: idx[i+1]:: " << idx[i + 1] << "\n";
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
			const Pos& p1 = H[p.hi][i0];
			std::cout << "hull rotate::\n";
			std::cout << "        cur::" << cur << "\n";
			std::cout << "          p::" << p << "\n";
			std::cout << "         p1::" << p1 << "\n";
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
			//std::cout << "outl swap:: i:: " << i << " j:: " << j << "\n";
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
		std::cout << "head[0]:: " << head[0] << "\n";
		vec = P[head[0]] - P[outl[K]];
		if (valid(cur, vec)) SQ.push(Event(vec, outl[K], head[0], HEAD));
		std::cout << "tail[0]:: " << tail[0] << "\n";
		vec = P[tail[0]] - P[outl[K]];
		if (valid(cur, vec)) SQ.push(Event(vec, outl[K], tail[0], TAIL));
		//std::cout << "loop start::\n";
		//std::cout << (t == BOT ? "BOT::\n" : "TOP::\n");
		//debug_pq(SQ, "SQ");
		while (SQ.size()) {
			Event ev = SQ.top();
			if (cur / ev.v > 0) break;
			SQ.pop();
			if (cur / ev.v < 0) continue;
			//std::cout << "ev.v:: " << ev.v << "\n";
			int i = ev.i, j = ev.j, typ = ev.t;
			std::cout << "i:: " << i << " j:: " << j << " type:: " << (typ == HEAD ? "HEAD" : "TAIL") << "\n";
			std::cout << "outl[K]:: " << outl[K] << "\n";
			Pos c, x = P[outl[K]];

			//이벤트가 머리에서 일어나는지 꼬리에서 일어나는지는 중요하지 않다.
			//결국 교환 과정은 앞뒤 모두에서 일어날 수 있기 때문.
			//후보점의 위치를 파악한 후 교환, 후보점군 갱신(머리, 꼬리 동시에)
			//튀어나온 점이 어디로 가야하는지를 파악해주어야 함
			//이게 좀 귀찮기는 함. 이 조건 분기를 생각해보자.

			std::cout << "x:: " << x << "\n";
			if (typ == HEAD) {
				if (outl[K] != ev.i || head[0] != ev.j) continue;
				c = P[head[0]];
			}
			else if (typ == TAIL) {
				if (outl[K] != ev.i || tail[0] != ev.j) continue;
				c = P[tail[0]];
			}
			EV.push_back(Event(cur, ord[x.pi][OUTL], -1, t));
			pop(x);
			std::cout << "pop::\n";
			int sz = H[c.hi].size();
			Pos pre = H[c.hi][(c.i - 1 + sz) % sz];
			Pos nxt = H[c.hi][(c.i + 1) % sz];
			std::cout << "pre, nxt::\n";
			std::cout << "c.pi  :: " << c.pi << "\n";
			std::cout << "c     :: " << c << "\n";
			std::cout << "pre.pi:: " << pre.pi << "\n";
			std::cout << "nxt.pi:: " << nxt.pi << "\n";
			std::cout << "sts[c.pi]:: " << sts[c.pi] << "\n";
			if (sts[c.pi] != -1 && sts[c.pi] & HEAD) {
			//if (sts[c.pi] & HEAD) {
				if (sts[nxt.pi] != OUTL) candidate_update(nxt, cur, HQ, head, HEAD, ord[c.pi][HEAD]);
				else candidate_update(nxt, cur, HQ, head, HEAD, ord[c.pi][HEAD], 1);
			}
			std::cout << "HEAD done\n";
			std::cout << "sts[c.pi]:: " << sts[c.pi] << "\n";
			if (sts[c.pi] != -1 && (sts[c.pi] & TAIL)) {
			//if (sts[c.pi] & TAIL) {
				if (sts[pre.pi] == OUTL && sts[nxt.pi] != OUTL && !(sts[nxt.pi] & TAIL)) {
					std::cout << "FUCK:: 1\n";
					candidate_update(nxt, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
				}
				else if (sts[pre.pi] != OUTL) {
					std::cout << "FUCK:: 2\n";
					std::cout << "              c:: " << c << "\n";
					std::cout << "ord[c.pi][TAIL]:: " << ord[c.pi][TAIL] << "\n";
					candidate_update(pre, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
				}
				else {
					std::cout << "FUCK:: 3\n";
					candidate_update(pre, cur, TQ, tail, TAIL, ord[c.pi][TAIL], 1);
				}
			}
			std::cout << "head, tail update::\n";

			for (int k = 0; k < hp; k++) {
				std::cout << "	" << (t == BOT ? "bot" : "top") << ".head[" << k << "]:: (" << P[head[k]].x << ", " << P[head[k]].y << ")\n";
			}
			for (int k = 0; k < tp; k++) {
				std::cout << "	" << (t == BOT ? "bot" : "top") << ".tail[" << k << "]:: (" << P[tail[k]].x << ", " << P[tail[k]].y << ")\n";
			}

			push(c);
			std::cout << "push::\n";
			vec = P[outl[K]] - P[outl[K - 1]];
			if (valid(cur, vec)) {
				OQ.push(Event(vec, outl[K - 1], outl[K]));
			}
			std::cout << "get vec::\n";
			sz = H[x.hi].size();
			pre = H[x.hi][(x.i - 1 + sz) % sz];
			nxt = H[x.hi][(x.i + 1) % sz];
			std::cout << "  x:: " << x << "\n";
			std::cout << "pre:: " << pre << " sts:: " << sts[pre.pi] << "\n";
			std::cout << "nxt:: " << nxt << " sts:: " << sts[nxt.pi] << "\n";
			if ((sts[x.pi] != -1 && !(sts[x.pi] & TAIL)) && cnt[x.hi] + 1 == sz) {
				std::cout << "TAIL 1:: FUCK::FUCK::FUCK::FUCK::FUCK::FUCK::\n";
				std::cout << "       sz:: " << sz << "\n";
				std::cout << "     x.hi:: " << x.hi << "\n";
				std::cout << "cnt[x.hi]:: " << cnt[x.hi] << "\n";
				candidate_update(x, cur, TQ, tail, TAIL, -1);
			}
			else if (sts[pre.pi] != -1 && (sts[pre.pi] & TAIL)) {
				std::cout << "TAIL 2:: FUCK::FUCK::FUCK::FUCK::FUCK::FUCK::\n";
				candidate_update(x, cur, TQ, tail, TAIL, ord[pre.pi][TAIL]);
			}
			vec = x - P[outl[K]];
			std::cout << "get x-vec::\n";
			std::cout << "cnt[x.hi]:: " << cnt[x.hi] << "\n";
			std::cout << "sz       :: " << sz << "\n";
			std::cout << "x        :: " << x << "\n";
			std::cout << "P[outl[K]]: " << P[outl[K]] << "\n";
			if ((sts[x.pi] != -1 && !(sts[x.pi] & HEAD)) && cnt[x.hi] + 1 == sz) {
				candidate_update(x, cur, HQ, head, HEAD, -1);
			}
			//else if ((sts[nxt.pi] | HEAD) && valid(cur, vec)) {
			else if (sts[nxt.pi] != -1 && (sts[nxt.pi] & HEAD)) {
				std::cout << "head prev update FUCK::FUCK::FUCK::FUCK::FUCK::FUCK::\n";
				candidate_update(x, cur, HQ, head, HEAD, ord[nxt.pi][HEAD]);
			}
			//else if (sts[pre.pi] | HEAD) candidate_update(x, cur, HQ, head, HEAD, ord[pre.pi][HEAD]);

			for (int k = 0; k < hp; k++) {
				std::cout << "	" << (t == BOT ? "bot" : "top") << ".head[" << k << "]:: (" << P[head[k]].x << ", " << P[head[k]].y << ")\n";
			}
			for (int k = 0; k < tp; k++) {
				std::cout << "	" << (t == BOT ? "bot" : "top") << ".tail[" << k << "]:: (" << P[tail[k]].x << ", " << P[tail[k]].y << ")\n";
			}

			//if (cnt[c.hi] == H[c.hi].size() - 1) {
			//	candidate_update(c, cur, HQ, head, HEAD, ord[c.pi][HEAD], 1);
			//	candidate_update(c, cur, TQ, tail, TAIL, ord[c.pi][TAIL], 1);
			//}
			//else if (cnt[c.hi] == H[c.hi].size() - 2) {
			//}
			//else {
			//	if (cnt[c.hi] <= H[c.hi].size() - 1 && sts[c.pi] | HEAD) candidate_update(nxt, cur, HQ, head, HEAD, ord[c.pi][HEAD]);
			//	if (cnt[c.hi] <= H[c.hi].size() - 1 && sts[c.pi] | TAIL) candidate_update(pre, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
			//}
			
			/*
			if (typ == HEAD) {
				if (outl[K] != ev.i || head[0] != ev.j) continue;
				c = P[head[0]];
			}
			else if (typ == TAIL) {
				if (outl[K] != ev.i || tail[0] != ev.j) continue;
				c = P[tail[0]];
			}
			int hi = c.hi;
			//std::cout << "hi:: " << hi << "\n";
			int sz = H[hi].size();
			const Pos& nxt = H[c.hi][(c.i + 1) % sz];
			//std::cout << "nxt:: " << nxt << "\n";
			//std::cout << "c.pi:: " << c.pi << "\n";
			if (sts[c.pi] | HEAD) candidate_update(nxt, cur, HQ, head, HEAD, ord[c.pi][HEAD]);
			else if (sts[c.pi] | TAIL) candidate_update(nxt, cur, TQ, tail, TAIL, ord[c.pi][TAIL]);
			//std::cout << "candi update done\n";
			pop(x);
			push(c);
			vec = P[outl[K]] - P[outl[K - 1]];
			if (valid(cur, vec)) {
				OQ.push(Event(vec, outl[K - 1], outl[K]));
			}
			//std::cout << "x pop and c push\n";
			if (cnt[x.hi] == H[x.hi].size()) {
				//std::cout << "cnt[x] == H[x].sz::\n";
				cnt[x.hi]--;
				candidate_update(x, cur, HQ, head, HEAD, -1);
				candidate_update(x, cur, TQ, tail, TAIL, -1);
			}
			else {
				//std::cout << "cnt[x] != H[x].sz!!::\n";
				cnt[x.hi]--;
				const Pos& pre = H[x.hi][(x.i - 1 + sz) % sz];
				if (sts[pre.pi] | TAIL) {
					candidate_update(x, cur, TQ, tail, TAIL, ord[pre.pi][TAIL]);
				}
				if (sts[pre.pi] | HEAD) {
					candidate_update(x, cur, HQ, head, HEAD, ord[pre.pi][HEAD]);
				}
			}
			*/
		}
		return;
	}
	bool jaw_rotate(Events& EV, const Pos& cur) {
		bool f = 0;
		//std::cout << "DEBUG:: cur:: " << cur << "\n";
		//debug_pq(HQ, "HQ");
		std::cout << (t == BOT ? "BOT " : "TOP ") << "candi head rotate::\n";
		candidate_jaw_rotate(cur, HQ, head, HEAD);
		//debug_pq(TQ, "TQ");
		std::cout << (t == BOT ? "BOT " : "TOP ") << "candi tail rotate::\n";
		candidate_jaw_rotate(cur, TQ, tail, TAIL);
		//debug_pq(OQ, "OQ");
		std::cout << (t == BOT ? "BOT " : "TOP ") << "outlier rotate::\n";
		outlier_jaw_rotate(cur, EV);
		//debug_pq(SQ, "SQ");
		std::cout << (t == BOT ? "BOT " : "TOP ") << "outlier candi swap::\n";
		outlier_candidate_swap(cur, EV);
		//debug_pq(JQ, "JQ");
		std::cout << (t == BOT ? "BOT " : "TOP ") << "hull rotate::\n";
		hull_jaw_rotate(cur);
		return f;
	}
};
struct Calipers {
	const Pos s = Pos(0, -1), e = Pos(0, 1);
	int N, K, h;
	Pos cur;
	Jaw bot = Jaw(BOT), top = Jaw(TOP);
	Calipers() { N = -1, K = -1, h = 0; }
	void init(int n = -1, int k = -1, Pos c = Pos(0, -1)) {
		N = n; K = k; cur = c;
		//std::sort(P.begin(), P.end(), cmpx_rvs);//굳이 역오름차순 정렬을 해야할 이유가?
		std::sort(P.begin(), P.end());//이 정렬 기준으로 정렬하면 초기에 무시해야하는 이벤트 기울기를 뛰어넘고 시작하는 게 가능함
		for (int i = 0; i < N; i++) P[i].pi = i;
		Polygon Q = P;
		Vbool F(N, 0);
		//std::cout << "init start::\n";
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
			for (int i = b; i <= sz; i++) U[k].push_back(H[k][i % sz]);
			for (const Pos& h : H[k]) F[h.pi] = 1;
			Q.clear();
			for (const Pos& t : tmp) if (!F[t.pi]) Q.push_back(t);
		}
		//std::cout << "hull done::\n";
		//std::sort(P.begin(), P.end(), cmpx_rvs);
		std::sort(P.begin(), P.end());
		for (int k = 0; k < h; k++) {
			for (const Pos& p : H[k]) P[p.pi].hi = p.hi, P[p.pi].i = p.i;
		}
		//std::cout << "numbering done::\n";
		Q.clear(); for (const Pos& p : P) if (F[p.pi]) Q.push_back(p);
		//for (int h_ = 0; h_ < h; h_++) {
		//	const Polygon& L = H[h_];
		//	std::cout << "Layer[" << h_ << "]::\n";
		//	for (const Pos& p : L) std::cout << "	" << p << "\n";
		//	std::cout << "Layer[" << h_ << "]::\n";
		//}
		//for (int h_ = 0; h_ < h; h_++) {
		//	const Polygon& LH = L[h_];
		//	std::cout << "LH[" << h_ << "]::\n";
		//	for (const Pos& p : LH) std::cout << "	" << p << "\n";
		//	std::cout << "LH[" << h_ << "]::\n";
		//	const Polygon& UH = U[h_];
		//	std::cout << "UH[" << h_ << "]::\n";
		//	for (const Pos& p : UH) std::cout << "	" << p << "\n";
		//	std::cout << "UH[" << h_ << "]::\n";
		//}
		bot.init(Q, N, K, h, s); top.init(Q, N, K, h, e);
		//std::cout << "bot, top done::\n";
		return;
	}
	bool jaw_rotate() {
		Polygon V;
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
		std::cout << "get vertical dist:::: b:: " << b << " t:: " << t << "\n";
		if (b < 0 || K < b || t < 0 || K < t) return INF;
		std::cout << "P[bot[outl[b]]:: " << P[bot.outl[b]] << " P[top.outl[t]]:: " << P[top.outl[t]] << "\n";
		return vertical_dist(P[bot.outl[b]], cur, P[top.outl[t]]);
	}
	bool rotate(ld& d) {
		static int CNT = 0;
		CNT++;
		int tq = sign(e / cur);
		//if (CNT == LIMIT) return 0;
		if (tq > 0 || (!tq && (e * cur) > 0)) return 0;
		//std::cout << "DEBUG:: TOP::\n";
		std::vector<Event> EV;
		std::cout << "\n\ncur.vec:: " << cur << "\n";
		//std::cout << "now rotate::\n";
		bot.jaw_rotate(EV, cur);
		//std::cout << "bot rotate::\n";
		top.jaw_rotate(EV, -cur);
		//std::cout << "top rotate::\n";
		//std::cout << "EV.size:: " << EV.size() << "\n";
		std::cout << "DEBUG:: cur:: (" << cur.x << ", " << cur.y << ")\n";
		std::cout << "DEBUG:: BOT::\n";
		for (int k = 0; k <= bot.K; k++) {
			std::cout << "	bot[" << k << "]:: (" << P[bot.outl[k]].x << ", " << P[bot.outl[k]].y << ") " << bot.sts[bot.outl[k]] << "\n";
		}
		for (int k = 0; k < bot.hp; k++) {
			std::cout << "	bot.head[" << k << "]:: (" << P[bot.head[k]].x << ", " << P[bot.head[k]].y << ") " << bot.sts[bot.head[k]] << "\n";
		}
		for (int k = 0; k < bot.tp; k++) {
			std::cout << "	bot.tail[" << k << "]:: (" << P[bot.tail[k]].x << ", " << P[bot.tail[k]].y << ") " << bot.sts[bot.tail[k]] << "\n";
		}
		for (int k = 0; k < bot.h; k++) {
			std::cout << "	bot.cnt[" << k << "]:: " << bot.cnt[k] << "\n";
		}
		std::cout << "DEBUG:: BOT::\n";
		std::cout << "DEBUG:: TOP::\n";
		for (int k = 0; k <= top.K; k++) {
			std::cout << "	top[" << k << "]:: (" << P[top.outl[k]].x << ", " << P[top.outl[k]].y << ") " << top.sts[top.outl[k]] << "\n";
		}
		for (int k = 0; k < top.hp; k++) {
			std::cout << "	top.head[" << k << "]:: (" << P[top.head[k]].x << ", " << P[top.head[k]].y << ") " << top.sts[top.head[k]] << "\n";
		}
		for (int k = 0; k < top.tp; k++) {
			std::cout << "	top.tail[" << k << "]:: (" << P[top.tail[k]].x << ", " << P[top.tail[k]].y << ") " << top.sts[top.tail[k]] << "\n";
		}
		for (int k = 0; k < top.h; k++) {
			std::cout << "	top.cnt[" << k << "]:: " << top.cnt[k] << "\n";
		}
		std::cout << "\n\n";
		for (const Event& ev : EV) {
			if (ev.t == BOT) d = std::min(d, dist(ev.i, K - ev.i));
			else d = std::min(d, dist(K - ev.i, ev.i));
		}
		//std::cout << "rotate rotate::\n";
		std::cout << "ret update:: " << d << "\n\n\n";
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
	std::cin >> N >> K;
	P.resize(N); for (Pos& p : P) std::cin >> p;
	if (N < 3) { std::cout << "0.000000000\n"; return; }
	if (N <= 3 && K == 1) { std::cout << "0.000000000\n"; return; }
	if (K == 0) { std::cout << rotating_calipers(P) << "\n"; return; }
	C.init(N, K, Pos(0, -1));
	//std::cout << "init done::\n";
	//int h = C.h;
	//for (int h_ = 0; h_ < h; h_++) {
	//	const Polygon& L = H[h_];
	//	std::cout << "Layer[" << h_ << "]::\n";
	//	for (const Pos& p : L) std::cout << "	" << p << "\n";
	//	std::cout << "Layer[" << h_ << "]::\n";
	//}
	ld ret = INF;
	//std::cout << "\n\n\n\n\n";
	//std::cout << "ret:: " << ret << "\n";
	//std::cout << "\n\n\n\n\n";
	while (C.rotate(ret)) {}
	std::cout << std::max(.0, ret * .5) << "\n";
	return;
}
int main() { solve(); return 0; }//boj26108


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

//볼록껍질 뒤쪽을 신경을 못 써줬다.
//구현하다 보니 참고했던 코드랑 비슷해져버리는 느낌이다.
//볼록껍질이 이벤트 라인과 걸쳐있는 점군들을 따로 한 쌍 씩 관리해서
//앞뒤쪽 갱신되는 후보들을 모두 관리할 수 있도록 한다.