#include <bits/stdc++.h>

using namespace std;

#define StarBurstStream ios_base::sync_with_stdio(false); cin.tie(0);
#define iter(a) a.begin(), a.end()
#define pb(a) emplace_back(a)
#define mp(a, b) make_pair(a, b)
#define ff first
#define ss second
#define topos(a) ((a) = (((a) % MOD + MOD) % MOD))
#define uni(a) a.resize(unique(iter(a)) - a.begin())
#define SZ(a) int(a.size())
#ifdef LOCAL
void debug() { cerr << "\n"; }
template<class T, class ... U> void debug(T a, U ... b) { cerr << a << " ", debug(b...); }
template<class T> void pary(T l, T r) {
	while (l != r) cerr << *l << " ", l++;
	cerr << "\n";
}
#else
#define debug(...) void()
#define pary(...) void()
#endif

typedef long long ll;
//typedef long double ld;

using pii = pair<int, int>;
using pll = pair<ll, ll>;

const ll MOD = 1000000007;
const ll MAX = 2147483647;

template<typename A, typename B>
ostream& operator<<(ostream& o, pair<A, B> p) {
	return o << '(' << p.ff << ',' << p.ss << ')';
}

ll ifloor(ll a, ll b) {
	if (b < 0) a *= -1, b *= -1;
	if (a < 0) return (a - b + 1) / b;
	else return a / b;
}

ll iceil(ll a, ll b) {
	if (b < 0) a *= -1, b *= -1;
	if (a > 0) return (a + b - 1) / b;
	else return a / b;
}

using Line = pair<pll, pll>;
using ld = long double;
using pdd = pair<ld, ld>;
#define X first
#define Y second
// ll eps = 1e-7;

pll operator+(pll a, pll b)
{
	return { a.X + b.X, a.Y + b.Y };
}
pll operator-(pll a, pll b)
{
	return { a.X - b.X, a.Y - b.Y };
}
pll operator-(pll a)
{
	return { -a.X, -a.Y };
}
pll operator*(ll i, pll v)
{
	return { i * v.X, i * v.Y };
}
pll operator*(pll v, ll i)
{
	return { i * v.X, i * v.Y };
}
pll operator/(pll v, ll i)
{
	return { v.X / i, v.Y / i };
}
ll dot(pll a, pll b)
{
	return a.X * b.X + a.Y * b.Y;
}
ll cross(pll a, pll b)
{
	return a.X * b.Y - a.Y * b.X;
}
ll abs2(pll v)
{
	return v.X * v.X + v.Y * v.Y;
};
ld abs(pll v)
{
	return sqrt(abs2(v));
};
int sgn(ll v)
{
	return v > 0 ? 1 : (v < 0 ? -1 : 0);
}
int ori(pll a, pll b, pll c)
{
	return sgn(cross(b - a, c - a));
}
bool collinearity(pll a, pll b, pll c)
{
	return ori(a, b, c) == 0;
}
bool btw(pll p, pll a, pll b)
{
	return collinearity(p, a, b) && sgn(dot(a - p, b - p)) <= 0;
}

bool seg_intersect(Line a, Line b) {
	pll p1, p2, p3, p4;
	tie(p1, p2) = a; tie(p3, p4) = b;
	if (btw(p1, p3, p4) || btw(p2, p3, p4) || btw(p3, p1, p2) || btw(p4, p1, p2))
		return true;
	return ori(p1, p2, p3) * ori(p1, p2, p4) < 0 &&
		ori(p3, p4, p1) * ori(p3, p4, p2) < 0;
}
pdd intersect(Line a, Line b) {
	pll p1, p2, p3, p4;
	tie(p1, p2) = a; tie(p3, p4) = b;
	ll a123 = cross(p2 - p1, p3 - p1);
	ll a124 = cross(p2 - p1, p4 - p1);
	pll tmp = p4 * a123 - p3 * a124;
	ll q = a123 - a124;
	//debug("intersect", p1, p2, p3, p4, pdd((ld)tmp.X / q, (ld)tmp.Y / q), tmp, q);
	return { (ld)tmp.X / q, (ld)tmp.Y / q };
}
pll perp(pll p1)
{
	return pll(-p1.Y, p1.X);
}
pll projection(pll p1, pll p2, pll p3)
{
	return p1 + (p2 - p1) * dot(p3 - p1, p2 - p1) / abs2(p2 - p1);
}

int is_neg(pll k, pll vec = pll(1, 0)) {
	return cross(vec, k) < 0 || (cross(vec, k) == 0 && dot(vec, k) < 0);
}
int cmp(pll a, pll b, pll vec = pll(1, 0), bool same = true) {
	int A = is_neg(a, vec), B = is_neg(b, vec);
	if (A != B)
		return A < B;
	if (sgn(cross(a, b)) == 0)
		return same ? abs2(a) < abs2(b) : -1;
	return sgn(cross(a, b)) > 0;
}

int n, K;

struct event {
	pll v;
	int ty = -1, lay, id;
	// 0: new point, 1: mn change, 2: layer rotate
	// 4: internal swap
};

ostream& operator<<(ostream& o, event e) {
	return o << '[' << e.v << ',' << e.ty << ',' << e.lay << ',' << e.id << ']';
}
bool operator<(event x, event y) {
	int t = cmp(x.v, y.v, pll(1, 0), false);
	//debug("comp1", x.v, y.v, t);
	if (t == -1) return x.ty > y.ty;
	return t;
}
bool operator>(event x, event y) {
	int t = cmp(x.v, y.v, pll(1, 0), false);
	//debug("comp1", x.v, y.v, t);
	if (t == -1) return x.ty < y.ty;
	return !t;
}

void print_pq(priority_queue<event, vector<event>, greater<>> pq) {
	while (!pq.empty()) {
		cerr << pq.top() << " ";
		pq.pop();
	}
	cerr << "\n";
}

ld ans = 1e18;

struct Sweep {
	vector<pll> pts;//점군
	vector<vector<pll>> poly;//볼록껍질 관리
	vector<event> knxt;//볼록껍질 선분 정렬 후 관리
	priority_queue<event, vector<event>, greater<>> pq;//매번 가장 기울기가 덜 변하는 점 쌍들을 관리하는 자료구조
	//이벤트 구조체를 순방향 정렬하는 pq, 2번째 인자는 원래 pq가 기본적으로 벡터로 구현되는데, 순서를 바꾸려면 c언어 특성상 앞 인자를 지정을 안 해주면 문제가 생겨서 넣었다고 함.
	vector<int> lp, rp, bot, top;//아직 모르겠음 ?????
	pll cur;//현재 기울기
	int clay = -1, cid = -1;//현재 가장 앞부분에 있는 점의 층과 볼록껍질 점 번호
	vector<pii> ord;//정렬된 순서, 껍질 번호 - 껍질 내 점 번호 순
	vector<vector<int>> pos;//이건 또 뭐야 ?????

	Sweep(vector<pll>& _p) : pts(_p), lp(K + 1, -1), rp(K + 1, -1), bot(K + 1, -1), top(K + 1, -1), cur(1, 0),
		ord(K + 1), pos(K + 1) {//모든 점군은 K + 1 개를 넘지 않도록 미리 정리되며, 이는 양쪽에서 최대 300개 점을 제외한 벨트의 폭을 바로 구하기 위함임
		poly.resize(K + 1);
		make_convex_layer();//볼록껍질 초기화
		debug("OwO");
		init();//?
	}

	bool good(pll p) {
		pll cp = poly[clay][cid];
		//ori 함수는 ccw 함수와 기능이 같음
		return ori(cp, cp + cur, p) > 0 || (ori(cp, cp + cur, p) == 0 && dot(cur, p - cp) <= 0);//현재 기울기보다 안쪽에 있거나 앞쪽에 있다
	}

	void upd_ans(int p, Sweep& another) {
		if (p < 0 || p > K) return;//점의 수가 넘어가면 벡터 조회에서 에러가 나며 계산을 할 수 있는 조건도 아님
		if (cross(cur, pll(0, 1)) == 0) return;//기울기가 마지막까지 도달한 상태
		pll a = getp(ord[p]);
		pll b = -another.getp(another.ord[K - p]);
		ld tans = abs(intersect(Line(a, a + cur), Line(b, b + pll(0, 1))).Y - b.Y);
		/*debug("upd_ans", cur, a, b, tans);
		for(pii i : ord) cerr << getp(i) << " ";
		cerr << "\n";
		for(pii i : another.ord) cerr << getp(i) << " ";
		cerr << "\n";*/
		ans = min(ans, tans / 2);
	}

	void make_convex_layer() {//볼록껍질 초기화
		debug("convex");
		vector<bool> use(n);
		poly.resize(K + 1);
		vector<int> tord;
		for (int i = 0; i < n; i++) tord.pb(i);
		sort(iter(tord), [&](int x, int y) { return pts[x] < pts[y]; });
		int r = n;
		for (int t = 0; t < K + 1; t++) {
			//if(t % 100 == 0) debug("test", t);
			vector<int> hull;
			for (int _ = 0; _ < 2; _++) {
				int sz = SZ(hull);
				for (int i : tord) {
					if (use[i]) continue;
					pll p = pts[i];
					while (SZ(hull) - sz >= 2 && ori(pts[hull.end()[-2]], pts[hull.back()], p) <= 0)
						hull.pop_back();
					hull.pb(i);
				}
				if (!hull.empty() && !(r == 1 && _ == 1)) hull.pop_back();
				reverse(iter(tord));
			}
			for (int i : hull) {
				use[i] = true;
				poly[t].pb(pts[i]);
				pary(iter(poly[t]));
				r--;
			}
			pos[t].resize(SZ(poly[t]), -1);
		}
		debug("convex ok");
	}

	pll get_nxt_vec() {//더 바깥쪽을 향하는 벡터를 선정하는 함수
		if (knxt.empty() && pq.empty()) return pll(0, 0);
		if (knxt.empty()) return pq.top().v;
		if (pq.empty()) return knxt.front().v;
		return min(knxt.front().v, pq.top().v, [&](pll x, pll y) { return cmp(x, y); });
	}

	void check_inter(int i) {
		if (i <= 0 || i > K) return;
		pll p1 = getp(ord[i]);//어떤 껍질 내의 어떤 점을 꺼냄
		pll p2 = getp(ord[i - 1]);//어떤 껍질 내의 어떤 점을 꺼냄
		//debug("check_inter", i, p1, p2);
		pll v = p2 - p1;//두 점 중 더 앞에 있는 점으로 향하는 벡터?
		if (!cmp(cur, v)) return;
		//debug("ok push", v);
		pq.push(event({ v, 4, ord[i].ff, ord[i].ss }));//2개의 이벤트로부터 만들어지는 기울기
	}

	void process_knxt(pll vec, Sweep& another) {
		if (knxt.empty() || cmp(vec, knxt.front().v, pll(1, 0), false) != -1) return;//다음 이벤트가 없거나 기울기가 존재한다면 return?
		/*debug("test", cur, "clay", clay, "cid", cid, "cp", poly[clay][cid]);
		for(int i = 0; i <= K; i++) cerr << getp(ord[i]) << " ";
		cerr << "\n";*/
		/*for(pll i : pts){
			if(good(i)) cerr << i << " ";
		}
		cerr << "\n";*/
		bool changed = false;
		for (auto [v, ty, lay, id] : knxt) {
			//debug("event", v, ty, lay, id);
			if (ty == 2) { // layer rotate
				int topnxt = getnxt(lay, top[lay]);//볼록껍질의 맨 위 점
				int botnxt = getnxt(lay, bot[lay]);//볼록껍질의 맨 아래 점
				pll vtop = poly[lay][top[lay]] - poly[lay][topnxt];//캘리퍼스 작동 준비
				pll vbot = poly[lay][botnxt] - poly[lay][bot[lay]];//캘리퍼스 작동 준비
				bool bot_rot = cmp(vbot, v, pll(1, 0), false) == -1;//pll(1,0)은 사분면을 반으로 나누기 위함. 내 정렬 방법이랑 비슷한데 구분선을 임의로 정해줄 수 있음
				bool top_rot = cmp(vtop, v, pll(1, 0), false) == -1;//두 벡터가 같은가?를 본다. 지금 캘리퍼스의 턱이 도착한 순간 보는 이 기울기가 일치한다면, 돌아간다.
				bool crot = false, allrot = false;
				if (top_rot && pii(lay, top[lay]) == pii(clay, cid)) crot = true;//이건 지금 보는 기울기와 뭔가를 비교하는 모습
				if (bot_rot && pii(lay, bot[lay]) == pii(clay, cid)) crot = true;
				if (top_rot && lp[lay] == rp[lay] && lp[lay] == top[lay]) allrot = true;//이건 뭘까?
				if (bot_rot && lp[lay] == rp[lay] && lp[lay] == bot[lay]) allrot = true;

				if (crot) {
					changed = true;
					if (top_rot && pii(lay, top[lay]) == pii(clay, cid)) {
						pos[clay][cid] = false;
						cid = getnxt(clay, cid);
						pos[clay][cid] = K;
						ord[K] = pii(clay, cid);
						upd_ans(K, another);
						check_inter(K);
					}
					else cid = getnxt(clay, cid);
				}

				if (allrot) {
					lp[lay] = getnxt(lay, lp[lay]);
					rp[lay] = getnxt(lay, rp[lay]);
				}

				if (bot_rot) bot[lay] = getnxt(lay, bot[lay]);
				if (top_rot) top[lay] = getnxt(lay, top[lay]);
			}
			else if (ty == 1) tie(clay, cid) = pii(lay, id); // mn change
			else if (!changed) { // ty == 0, new point
				if (top[clay] == cid) lp[clay] = rp[clay] = -1;
				else {
					if (cid == lp[clay]) lp[clay] = getlst(clay, lp[clay]);
					if (cid == rp[clay]) rp[clay] = getnxt(clay, rp[clay]);
				}
				pos[clay][cid] = -1;
				pos[lay][id] = K;
				ord[K] = pii(lay, id);
				upd_ans(K, another);
				check_inter(K);
				clay = lay; cid = id;
				if (top[lay] == id) lp[lay] = rp[lay] = id;
				else {
					if (getnxt(lay, lp[lay]) == id) lp[lay] = id;
					if (getlst(lay, rp[lay]) == id) rp[lay] = id;
				}
			}
		}
		find_knxt();
	}

	void process_pq(pll vec, Sweep& another) {
		if (pq.empty()) return;
		//debug("process_pq", pq.top().v);
		//print_pq(pq);
		while (!pq.empty() && cmp(pq.top().v, vec, pll(1, 0), false) == -1) {
			auto [v, ty, lay, id] = pq.top();
			pq.pop();
			int p = pos[lay][id];
			if (p <= 0) continue;
			pll p1 = poly[lay][id];
			auto [lay2, id2] = ord[p - 1];
			pll p2 = getp(ord[p - 1]);
			if (cmp(p2 - p1, v, pll(1, 0), false) != -1) continue;
			//debug("event", v, ty, lay, id);
			swap(ord[p - 1], ord[p]);
			swap(pos[lay][id], pos[lay2][id2]);
			upd_ans(p, another);
			upd_ans(p - 1, another);
			check_inter(p + 1);
			check_inter(p);
			check_inter(p - 1);
		}
	}

	void process(pll v, Sweep& another) {
		if (cmp(get_nxt_vec(), v, pll(1, 0), false) != -1) return;
		cur = v;
		process_knxt(v, another);
		process_pq(v, another);

		/*debug("process done");
		debug("test", cur, "clay", clay, "cid", cid, "cp", poly[clay][cid]);
		for(int i = 0; i <= K; i++) cerr << getp(ord[i]) << " ";
		cerr << "\n";
		for(pll i : pts){
			if(good(i)) cerr << i << " ";
		}
		cerr << "\n";*/
	}

	int getnxt(int lay, int id) {//볼록껍질의 다음 점을 찾는 함수
		return (id + 1) % SZ(poly[lay]);
	}
	int getlst(int lay, int id) {//볼록껍질의 이전 점을 찾는 함수
		return (id - 1 + SZ(poly[lay])) % SZ(poly[lay]);
	}
	pll getp(pii x) {
		auto [lay, id] = x;
		return poly[lay][id];//pii x 에는 볼록껍질의 번호와 그 껍질 내 점 번호가 들어있음
	}

	void init() {
		//debug("init");
		vector<pii> tmp;
		for (int i = 0; i <= K; i++)
			for (int j = 0; j < SZ(poly[i]); j++)
				tmp.pb(pii(i, j));

		sort(iter(tmp), [&](pii x, pii y) {
			pll px = getp(x), py = getp(y);
			return px.Y != py.Y ? px.Y > py.Y : px.X < py.X;
			});
		for (int i = 0; i <= K; i++) {
			ord[i] = tmp[i];
			pos[tmp[i].ff][tmp[i].ss] = i;
		}
		for (int i = 1; i <= K; i++) check_inter(i);

		tie(clay, cid) = tmp[K];
		auto good = [&](pll p) {
			pll cp = poly[clay][cid];
			return ori(cp, cp + cur, p) > 0 || (ori(cp, cp + cur, p) == 0 && dot(cur, p - cp) <= 0);
			};

		for (auto [lay, id] : tmp) {
			if (bot[lay] == -1) top[lay] = id;
			bot[lay] = id;
		}
		for (auto [lay, id] : tmp) {
			if (!good(poly[lay][id])) continue;
			//debug("owo", lay, id, poly[lay][id], ori(poly[lay][bot[lay]], poly[lay][top[lay]], poly[lay][id]));
			int t = ori(poly[lay][bot[lay]], poly[lay][top[lay]], poly[lay][id]);
			if (t >= 0) lp[lay] = id;
			if (t <= 0) rp[lay] = id;
		}
		debug("qAQ");

		find_knxt();
	}

	void find_knxt() {
		knxt.clear();

		vector<event> ev;
		pll cp = getp(pii(clay, cid));

		//debug("test", cur, "clay", clay, "cid", cid, "cp", cp);
		//debug("find_knxt");
		//for(int i = 0; i <= K; i++) debug("lay", i, bot[i], top[i], lp[i], rp[i]);
		/*for(pll i : pts){
			if(good(i)) cerr << i << " ";
		}
		cerr << "\n";*/

		auto mnchange = [&](int lay, int id) {
			if (pii(lay, id) == pii(clay, cid)) return;
			pll p = poly[lay][id];
			if (!good(p)) return;
			ev.pb(event({ p - cp, 1, lay, id }));
			};
		auto newpoint = [&](int lay, int id) {
			pll p = poly[lay][id];
			if (good(p)) return;
			ev.pb(event({ cp - p, 0, lay, id }));
			};
		auto rotate_layer = [&](int lay) {
			if (SZ(poly[lay]) < 2) return;
			int topnxt = getnxt(lay, top[lay]);
			int botnxt = getnxt(lay, bot[lay]);
			pll vtop = poly[lay][top[lay]] - poly[lay][topnxt];
			pll vbot = poly[lay][botnxt] - poly[lay][bot[lay]];
			pll v = cmp(vtop, vbot, cur) ? vtop : vbot;
			ev.pb(event({ v, 2, lay, -1 }));
			};

		for (int i = 0; i <= K; i++) {
			if (poly[i].empty()) continue;
			rotate_layer(i);
			if (lp[i] == -1) newpoint(i, top[i]);
			else {
				newpoint(i, getlst(i, rp[i]));
				if (getlst(i, rp[i]) != getnxt(i, lp[i]))
					newpoint(i, getnxt(i, lp[i]));
				mnchange(i, lp[i]);
				if (lp[i] != rp[i])
					mnchange(i, rp[i]);
			}
		}

		if (ev.empty()) return;
		//pary(iter(ev));
		pll nxt = ev[0].v;
		for (auto e : ev) {
			pll v = e.v;
			if (cmp(v, nxt, cur)) nxt = v;
		}
		if (is_neg(cur) && !is_neg(nxt)) return;

		for (auto [v, ty, lay, id] : ev) {
			if (cmp(v, nxt, pll(1, 0), false) == -1)
				knxt.pb(event({ v, ty, lay, id }));
		}
	}
};

int main() {
	StarBurstStream;

	cin >> n >> K;
	vector<pll> pts(n), pts2(n);
	for (int i = 0; i < n; i++) {
		cin >> pts[i].X >> pts[i].Y;
		pts2[i] = -pts[i];
	}
	debug("OK");

	//debug("init up");
	Sweep sup(pts);
	debug("OK");
	//debug("init down");
	Sweep sdown(pts2);
	debug("OK");
	for (int i = 0; i <= K; i++)
		sup.upd_ans(i, sdown);

	debug("OK");
	while (true) {
		pll vup = sup.get_nxt_vec();
		pll vdown = sdown.get_nxt_vec();
		pll v;
		if (vup == pll(0, 0)) v = vdown;
		else if (vdown == pll(0, 0)) v = vup;
		else v = min(vdown, vup, [&](pll x, pll y) { return cmp(x, y); });
		if (v == pll(0, 0)) break;
		//debug("do", v);

		//debug("do up", v);
		sup.process(v, sdown);
		//debug("do down", v);
		sdown.process(v, sup);
	}

	cout << fixed << setprecision(20) << ans << "\n";

	return 0;
}