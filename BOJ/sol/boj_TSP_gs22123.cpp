//THIS IS NOT MY CODE.
//AUTHOR: gs22123

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <random>
#include <array>
#include <tuple>
#include <multiset>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;
typedef pair<ll,ll> pll;
typedef pair<double,double> pdd;
#define fastio cin.tie(0)->sync_with_stdio(0); cout.tie(0);
#define all(x) x.begin(),x.end()
#define compress(x) x.erase(unique(all(x)),x.end())
#define ff first
#define ss second
#define MAX 500010
#define SIZE 100010

# pragma GCC optimize ("O3")
//# pragma GCC optimize ("Ofast")
# pragma GCC optimize ("unroll-loops")

random_device rd;
mt19937 gen(time(NULL));

namespace FASTIO{
	inline int readChar();
	template<class T = int> inline T readInt();
	template<class T> inline void writeInt(T x, char end = 0);
	inline void writeChar(int x);
	inline void writeWord(const char *s);
	static const int buf_size = 1 << 18;
	inline int getChar(){
	    #ifndef LOCAL
	    static char buf[buf_size];
	    static int len = 0, pos = 0;
	    if(pos == len) pos = 0, len = fread(buf, 1, buf_size, stdin);
	    if(pos == len) return -1;
	    return buf[pos++];
	    #endif
	}
	inline int readChar(){
	    #ifndef LOCAL
	    int c = getChar();
	    while(c <= 32) c = getChar();
	    return c;
	    #else
	    char c; cin >> c; return c;
	    #endif
	}
	template <class T>
	inline T readInt(){
	    #ifndef LOCAL
	    int s = 1, c = readChar();
	    T x = 0;
	    if(c == '-') s = -1, c = getChar();
	    while('0' <= c && c <= '9') x = x * 10 + c - '0', c = getChar();
	    return s == 1 ? x : -x;
	    #else
	    T x; cin >> x; return x;
	    #endif
	}
	static int write_pos = 0;
	static char write_buf[buf_size];
	inline void writeChar(int x){
	    if(write_pos == buf_size) fwrite(write_buf, 1, buf_size, stdout), write_pos = 0;
	    write_buf[write_pos++] = x;
	}
	template <class T>
	inline void writeInt(T x, char end){
	    if(x < 0) writeChar('-'), x = -x;
	    char s[24]; int n = 0;
	    while(x || !n) s[n++] = '0' + x % 10, x /= 10;
	    while(n--) writeChar(s[n]);
	    if(end) writeChar(end);
	}
	inline void writeWord(const char *s){
	    while(*s) writeChar(*s++);
	}
	struct Flusher{
	    ~Flusher(){ if(write_pos) fwrite(write_buf, 1, write_pos, stdout), write_pos = 0; }
	}flusher;
}

typedef long double lld;
const lld EPS = 1e-12, INF = 1e100;

struct P{
	double x,y;
	P(lld x = 0, lld y = 0):x(x),y(y){}

	P operator + (const P&p) const {
		return {x+p.x,y+p.y};
	}
	P operator - (const P&p) const {
		return {x-p.x,y-p.y};
	}
	P operator * (const lld&d) const {
		return {x*d,y*d};
	}
	P operator / (const ll&d) const {
		return {x/d,y/d};
	}
	bool operator != (const P&p) const {
		return x!=p.x||y!=p.y;
	}

	P rot() const {
		return P(-y,x);
	}

	lld dot(const P&p) const {
		return x*p.x + y*p.y;
	}
	lld cross(const P&p) const {
		return x*p.y-y*p.x;
	}

	lld len() const {
		return hypot(x,y);
	}

	bool operator < (const P&p) const {
		return make_pair(x,y)<make_pair(p.x,p.y);
	}
};

namespace Delaunator{
	bool collinear(P a, P b){
		return abs(a.cross(b)) < EPS;
	}

	P lineline(P a, P b, P c, P d){
		return a+(b-a)*((c-a).cross(d-c)/(b-a).cross(d-c));
	}

	P circumcenter(P a, P b, P c){
		b = (a+b)*0.5;
		c = (a+c)*0.5;
		return lineline(b,b+(b-a).rot(),c,c+(c-a).rot());
	}

	lld sweepx;

	struct arc{
		mutable P p,q;
		mutable int id = 0, i;
		arc(P p, P q, int i): p(p),q(q),i(i){}

		lld gety(lld x) const {
			if(q.y==INF) return INF;
			x += EPS;
			P med = (p+q)*0.5;
			P dir = (p-med).rot();
			lld D = (x-p.x)*(x-q.x);
			return med.y+((med.x-x)*dir.x+sqrtl(D)*dir.len())/dir.y;
		}
		bool operator<(const lld &y) const {
			return gety(sweepx) < y;
		}
		bool operator<(const arc &o) const {
			return gety(sweepx)<o.gety(sweepx);
		}
	};

	using beach = multiset<arc,less<> >;

	struct event{
		lld x;
		int id;
		beach::iterator it;
		event(lld x, int id, beach::iterator it):x(x),id(id),it(it){}

		bool operator<(const event&e) const {
			return x>e.x;
		}
	};

	struct fortune{
		beach line;
		vector<pair<P, int> > v;
		priority_queue<event> Q;
		vector<pii> edges;
		vector<bool> valid;
		int n, ti;
		fortune(vector<P> p){
			n = p.size();
			v.resize(n);
			for(int i=0;i<n;i++) v[i] = {p[i],i};
			sort(all(v));
		}

		void upd(beach::iterator it){
			if(it->i==-1) return;
			valid[-it->id] = false;
			auto a = prev(it);
			if(collinear(it->q-it->p,a->p-it->p)) return;
			it->id = --ti;
			valid.push_back(true);
			P c = circumcenter(it->p,it->q,a->p);
			lld x = c.x+(c-it->p).len();

			if(x>sweepx-EPS&&a->gety(x)+EPS>it->gety(x)){
				Q.push(event(x,it->id,it));
			}
		}

		void add_edge(int i, int j){
			if(i==-1||j==-1) return;
			edges.push_back({v[i].ss,v[j].ss});
		}

		void add(int i){
			P p = v[i].ff;
			auto c = line.lower_bound(p.y);
			auto b = line.insert(c,arc(p,c->p,i));
			auto a = line.insert(b,arc(c->p,p,c->i));
			add_edge(i,c->i);
			upd(a); upd(b); upd(c);
		}

		void remove(beach::iterator it){
			auto a = prev(it);
			auto b = next(it);
			line.erase(it);
			a->q = b->p;
			add_edge(a->i,b->i);
			upd(a); upd(b);
		}

		void solve(lld X = 1e9){
			X *= 3;
			line.insert(arc(P(-X,-X),P(-X,X),-1));
			line.insert(arc(P(-X,X),P(INF,INF),-1));

			for(int i=0;i<n;i++){
				Q.push(event(v[i].ff.x,i,line.end()));
			}
			ti = 0;
			valid.assign(1,false);
			while(!Q.empty()){
				event e = Q.top(); Q.pop();
				sweepx = e.x;
				if(e.id>=0){
					add(e.id);
				}
				else if(valid[-e.id]){
					remove(e.it);
				}
			}
		}
	};
}

ll n,k,centr[1000],visc[10000];
P pnt[10000],cent[1000];
vector<ll> clust[1000],g[1000];
vector<vector<ll> > ans;
vector<P> clustp[1000];
set<ll> ins[1000];

lld dist(P a, P b){
	return 1.0*hypot(a.x-b.x,a.y-b.y);
}

lld sqdist(P a, P b){
	return (a-b).x*(a-b).x+(a-b).y*(a-b).y;
}

double ddist(P a, P b){
	return 1.0*hypot(a.x-b.x,a.y-b.y);
}

lld dist(ll a, ll b){
	return dist(pnt[a],pnt[b]);
}

lld leng(vector<ll> v){
	ll sz = v.size();
	lld d = 0;
	for(int i=0;i<sz;i++){
		d += dist(v[i],v[(i+1)%sz]);
	}
	return d;
}

void getcents(){
	uniform_int_distribution<int> uid(0, (int)2e9);
	ll fir = uid(gen)%k+1;
	visc[fir] = 1;
	cent[1] = pnt[fir]; centr[1] = fir;
	for(int i=2;i<=k;i++){
		double d = 0;
		ll id = 1;
		for(int j=1;j<=n;j++){
			if(visc[j]) continue;
			lld dd = 1e30;
			for(int t=1;t<i;t++){
				dd = min(dd,dist(cent[t],pnt[j]));
			}
			if(dd>d){
				d = dd;
				id = j;
			}
		}
		cent[i] = pnt[id];
		centr[i] = id;
		visc[id] = 1;
	}
}

void kmeans(){
	ll nw = 1, cnt = 0;
	while(nw&&cnt<10000){
		P sum[1000] = {};
		for(int i=1;i<=k;i++){
			clust[i].clear();
			//clust[i].push_back(centr[i]);
		}
		for(int i=1;i<=n;i++){
			//if(visc[i]) continue;
			lld d = 1e30;
			ll idx = 1;
			for(int j=1;j<=k;j++){
				lld nd = sqdist(pnt[i],cent[j]);
				//if(i==centr[j]) cout << j << " " << nd << " " << d << " ";
				if(nd+EPS<=d){
					d = nd;
					idx = j;
				}
			}
			clust[idx].push_back(i);
			sum[idx] = sum[idx] + pnt[i];
		}
		nw = 0;
		for(int i=1;i<=k;i++){
			P nc = sum[i]/(ll)clust[i].size();
			//cout << clust[i].size() << "\n";
			//cout << nc.x << "a" << nc.y << " " << cent[i].x << " " << cent[i].y << "\n";
			if(nc!=cent[i]) nw++;
			cent[i] = nc;
		}
		cnt++;
		//cout << nw << "\n";
	}
}

vector<ll> initp(ll sz){
	uniform_int_distribution<int> uid(0, (int)2e9);
	ll vis[100010]={};
	vector<ll> in;
	for(int i=0;i<sz;i++){
		ll rd = uid(gen)%sz;
		while(vis[rd]) rd = uid(gen)%sz;
		in.push_back(rd);
		vis[rd] = 1;
	}
	return in;
}

vector<ll> nn(ll n, vector<P> p){
    uniform_int_distribution<int> uid(0, (int)2e9);
	ll vis[10000]={}, cnt = 1, prv = uid(gen)%n;
	vector<ll> path; path.push_back(prv); vis[prv] = 1;
	while(cnt<n){
		double d = DBL_MAX;
		ll id = 0;
		for(int i=0;i<n;i++){
			if(vis[i]) continue;
			if(d>ddist(p[prv],p[i])){
				d = ddist(p[prv],p[i]);
				id = i;
			}
		}
		vis[id] = 1;
		path.push_back(id);
		prv = id;
		cnt++;
	}
	return path;
}

lld rever(vector<ll> &path, ll i, ll j, ll k, vector<P>& p){
	ll sz = path.size();
	i %= sz; j %= sz; k %= sz;
	P i1 = p[path[i]], j1 = p[path[j]], k1 = p[path[k]];
	P i2 = p[path[(i+1)%sz]],j2 = p[path[(j+1)%sz]],k2 = p[path[(k+1)%sz]];
	vector<pair<double,ll> > tmp;
	tmp.push_back({ddist(i1,i2)+ddist(j1,j2)+ddist(k1,k2),1});
	tmp.push_back({ddist(i1,i2)+ddist(j1,k1)+ddist(j2,k2),2});
	tmp.push_back({ddist(i1,j1)+ddist(i2,j2)+ddist(k1,k2),3});
	tmp.push_back({ddist(i1,k1)+ddist(i2,k2)+ddist(j1,j2),4});
	tmp.push_back({ddist(i1,j1)+ddist(i2,k1)+ddist(j2,k2),5});
	tmp.push_back({ddist(i1,k1)+ddist(j2,i2)+ddist(j1,k2),6});
	tmp.push_back({ddist(i1,j2)+ddist(k1,j1)+ddist(i2,k2),7});
	tmp.push_back({ddist(i1,j2)+ddist(k1,i2)+ddist(j1,k2),8});
	sort(all(tmp));
	if(tmp[0].ss==1) return 10;
	ll D = tmp[0].ss;
	ll dif = tmp[0].ff - (ddist(i1,i2)+ddist(j1,j2)+ddist(k1,k2));

	if(D==2){
		for(int t=j+1;t<=(j+k+1)/2;t++){
			swap(path[t%sz],path[(j+k+1-t)%sz]);
		}
	}
	else if(D==3){
		for(int t=i+1;t<=(i+j+1)/2;t++){
			swap(path[t%sz],path[(i+j+1-t)%sz]);
		}
	}
	else if(D==4){
		for(int t=i+1;t<=(i+k+1)/2;t++){
			swap(path[t%sz],path[(i+k+1-t)%sz]);
		}
	}
	else if(D==5){
		for(int t=i+1;t<=(i+j+1)/2;t++){
			swap(path[t%sz],path[(i+j+1-t)%sz]);
		}
		for(int t=j+1;t<=(j+k+1)/2;t++){
			swap(path[t%sz],path[(j+k+1-t)%sz]);
		}
	}
	else if(D==6){
		for(int t=i+1;t<=(i+j+1)/2;t++){
			swap(path[t%sz],path[(i+j+1-t)%sz]);
		}
		for(int t=i+1;t<=(i+k+1)/2;t++){
			swap(path[t%sz],path[(i+k+1-t)%sz]);
		}
	}
	else if(D==7){
		for(int t=j+1;t<=(j+k+1)/2;t++){
			swap(path[t%sz],path[(j+k+1-t)%sz]);
		}
		for(int t=i+1;t<=(i+k+1)/2;t++){
			swap(path[t%sz],path[(i+k+1-t)%sz]);
		}
	}
	else if(D==8){
		for(int t=i+1;t<=(i+j+1)/2;t++){
			swap(path[t%sz],path[(i+j+1-t)%sz]);
		}
		for(int t=j+1;t<=(j+k+1)/2;t++){
			swap(path[t%sz],path[(j+k+1-t)%sz]);
		}
		for(int t=i+1;t<=(i+k+1)/2;t++){
			swap(path[t%sz],path[(i+k+1-t)%sz]);
		}
	}

	return dif;
}

vector<ll> opt3(vector<ll>& path, vector<P>& pn, ll thr){
	ll cnt = 0, n = path.size();
	while(cnt<thr){
		double del = 0;
		for(int i=0;i<n;i++){
			for(int j=i+2;j<n;j++){
				for(int k=j+2;k<n;k++){
					del += rever(path,i,j,k,pn);
					cnt++;
					if(cnt>thr) return path;
				}
			}
		}
		if(del+EPS>=0) break;
	}
	return path;
}

vector<ll> opt2(vector<ll> path, vector<P> pn, ll thr){
	ll cnt = 0,f=1,n=path.size();
	lld d = 0;
	for(int i=0;i<n;i++) d += dist(pn[path[i]],pn[path[(i+1)%n]]);
	while(f&&cnt<thr){
		f=0; cnt++;
		for(int i=0;i<=n-2;i++){
			for(int j=i+1;j<=n-1;j++){
				lld nl = -dist(pn[path[i]],pn[path[(i+1)%n]])-dist(pn[path[j]],pn[path[(j+1)%n]])
							+dist(pn[path[i]],pn[path[j]])+dist(pn[path[(i+1)%n]],pn[path[(j+1)%n]]);

				if(nl<0){
					//cout << nl << "\n";
					reverse(path.begin()+i+1,path.begin()+j+1);
					d += nl;
					f=1;
				}
			}
		}
	}
	return path;
}

struct pp{
	ll x,y,grp;
};

ll visp[1000][1000];
pp ansp[100000];
vector<double> G_dist(1000);
vector<ll> print_(1000);

ll pre[10000],suc[10000],grp[10000],siz[1000];
double tdist[1000];
set<pair<double,ll>,greater<pair<double,ll>> > st;

void nswap(ll i, ll j){
	if(siz[grp[i]]==1&&siz[grp[j]]==1) return;
	ll pi = pre[i], si = suc[i];
	ll pj = pre[j], sj = suc[j];
	if(siz[grp[i]]==1){
		pre[j] = j; suc[j] = j;
		pre[i] = pj; suc[i] = sj;
		suc[pj] = i; pre[sj] = i;
	}
	else if(siz[grp[j]]==1){
		pre[i] = i; suc[i] = i;
		pre[j] = pi; suc[j] = si;
		suc[pi] = j; pre[si] = j;
	}
	else{
		pre[i] = pj; suc[i] = sj;
		suc[pi] = j; pre[si] = j;
		pre[j] = pi; suc[j] = si;
		suc[pj] = i; pre[sj] = i;
	}
	swap(grp[i],grp[j]);
}

void insertp(ll i, ll j){
	if(siz[grp[j]]==1) return;
	siz[grp[j]]--; siz[grp[i]]++;
	ll pi = pre[i];
	ll pj = pre[j], sj = suc[j];
	suc[pj] = sj; pre[sj] = pj;
	suc[pi] = j; pre[i] = j;
	suc[j] = i; pre[j] = pi;
	grp[j] = grp[i];
}

void inserts(ll i, ll j){
	if(siz[grp[j]]==1) return;
	siz[grp[j]]--; siz[grp[i]]++;
	ll si = suc[i];
	ll pj = pre[j], sj = suc[j];
	suc[pj] = sj; pre[sj] = pj;
	suc[i] = j; pre[si] = j;
	suc[j] = si; pre[j] = i;
	grp[j] = grp[i];
}

void m2opt(ll thr){
	ll cnt = 0;
	for(ll i=1;i<=k;i++) for(auto j:g[i]){
		if(siz[i]==1&&siz[j]==1) continue;
		for(auto p1:ins[i]) for(auto p2:ins[j]){
			if(grp[p1]==grp[p2]) continue;
			double ch1=0,ch2=0;

			//swap
			if(siz[grp[p1]]==1){
				ch1 = 0;
				ch2 = dist(p1,suc[p2])+dist(p1,pre[p2])-dist(p2,suc[p2])-dist(p2,pre[p2]);
			}
			else if(siz[grp[p2]]==1){
				ch1 = dist(p2,suc[p1])+dist(p2,pre[p1])-dist(p1,suc[p1])-dist(p1,pre[p1]);
				ch2 = 0;
			}
			else{
				ch1 = dist(p2,suc[p1])+dist(p2,pre[p1])-dist(p1,suc[p1])-dist(p1,pre[p1]);
				ch2 = dist(p1,suc[p2])+dist(p1,pre[p2])-dist(p2,suc[p2])-dist(p2,pre[p2]);
			}
			if(max(tdist[grp[p1]],tdist[grp[p2]])>max(tdist[grp[p1]]+ch1,tdist[grp[p2]]+ch2)){
				tdist[grp[p1]] += ch1;
				tdist[grp[p2]] += ch2;
				nswap(p1,p2);
			}

			if(siz[grp[p2]]>1){
				//insertp1
				ch1 = dist(pre[p1],p2)+dist(p1,p2)-dist(pre[p1],p1);
				ch2 = dist(suc[p2],pre[p2])-dist(p2,pre[p2])-dist(p2,suc[p2]);
				if(max(tdist[grp[p1]],tdist[grp[p2]])>max(tdist[grp[p1]]+ch1,tdist[grp[p2]]+ch2)){
					tdist[grp[p1]] += ch1;
					tdist[grp[p2]] += ch2;
					insertp(p1,p2);
				}
			}

			if(siz[grp[p2]]>1){
				//inserts1
				ch1 = dist(suc[p1],p2)+dist(p1,p2)-dist(suc[p1],p1);
				ch2 = dist(suc[p2],pre[p2])-dist(p2,pre[p2])-dist(p2,suc[p2]);
				if(max(tdist[grp[p1]],tdist[grp[p2]])>max(tdist[grp[p1]]+ch1,tdist[grp[p2]]+ch2)){
					tdist[grp[p1]] += ch1;
					tdist[grp[p2]] += ch2;
					inserts(p1,p2);
				}
			}
			cnt++;
			if(cnt>thr) return;
		}
	}
}

int main(){
	fastio;
	ll k1 = 14, k2 = 10;
	ll tx = 1 + 814000 / k1, ty = 1 + 814000 / k2;
	uniform_int_distribution<int> uid(0, (int)2e9);
	n = FASTIO::readInt();
	k = FASTIO::readInt();
	ans.resize(k);
	vector<ll> v2[200][200] = {};
	for(int i=1;i<=n;i++){
		pnt[i].x = FASTIO::readInt();
		pnt[i].y = FASTIO::readInt();
		v2[(ll)pnt[i].x/tx][(ll)pnt[i].y/ty].push_back(i);
	}
	// ll vis[100010]={};
	// for(int i=1;i<=k;i++){
	// 	ll rd = uid(gen)%n+1;
	// 	while(vis[rd]) rd = uid(gen)%n+1;
	// 	cent[i] = pnt[rd];
	// 	vis[rd] = 1;
	// 	//cout << rd << "\n";
	// }
	// getcents();
	// kmeans();
	ll gr = 1;
	for(int i=0;i<k1;i++){
		for(int j=0;j<k2;j++){
			for(auto g:v2[i][j]) clust[gr].push_back(g);
			gr++;
		}
	}
	for(int i=1;i<=k;i++){
		for(auto j:clust[i]) clustp[i].push_back(pnt[j]);

		/*
		ACO::init(clustp[i]);
		vector<ll> path = ACO::aco(1,20);
		*/

		//DLAStsp::init(clustp[i]);
		ll sz = clustp[i].size();
		vector<ll> path = nn(sz,clustp[i]);
        //vector<ll> path = initp(sz);
		//path = DLAStsp::DLAS(path,100000);
		
		//path = opt2(path,clustp[i],100000);
		path = opt3(path,clustp[i],800000);
		//cout << path.size() << " ";
		ll path_size = path.size();
		siz[i] = path_size;
		for(ll j=0;j<path_size;j++){
			ans[i-1].push_back(clust[i][j]);

			ll curr=clust[i][path[j]];
			ll left=clust[i][path[(j-1+path_size)%path_size]];
			ll right=clust[i][path[(j+1+path_size)%path_size]];
			ansp[curr]={left,right,i};
			G_dist[i]+=dist(pnt[left],pnt[curr]);

			pre[curr] = left; suc[curr] = right; grp[curr] = i;
			ins[i].insert(curr);
        }
        tdist[i] = G_dist[i];

        // cout << G_dist[i] <<" ";
		// cout << path.size() << " ";
		// for(auto j:path) cout << clust[i][j] << " ";

		// cout << "\n";
	}
	//cout << "\n";

	vector<P> vv; vv.push_back(P(0,0));
	for(int i=1;i<=n;i++) vv.push_back(pnt[i]);
	// DLASmtsp::init(vv);
	// ans = DLASmtsp::DLAS(ans,20000);

	vector<P> tmp,tmp2;
	for(int i=1;i<=k;i++){
		P p = cent[i];
		tmp.push_back(P(p.x*cos(1)-p.y*sin(1),p.x*sin(1)+p.y*cos(1)));
	}
	Delaunator::fortune tr(tmp);
	tr.solve();
	vector<pii> dt = tr.edges;

	for(auto i:dt){
		g[i.ff+1].push_back(i.ss+1);
		g[i.ss+1].push_back(i.ff+1);
	}


	ll thr = 5500000;
	m2opt(thr);

	vector<ll> final;
	ll chk[10000] = {};
	for(int i=1;i<=n;i++){
		double d = 0;
		final.clear();
		if(chk[i]) continue;
		chk[i] = 1;
		final.push_back(i);
		ll pt = suc[i];
		d += dist(i,pt);
		while(!chk[pt]){
			final.push_back(pt);
			chk[pt] = 1;
			d += dist(pt,suc[pt]);
			pt = suc[pt];
		}

		//final = opt2(final,vv,200000);
        final = opt3(final,vv,25000);

		//cout << leng(final) << " ";
		//cout << tdist[grp[i]] << " ";
		//cout << tdist[grp[i]] << " " << d << "\n";
		FASTIO::writeInt(final.size(),' ');
		//cout << final.size() << " ";
		for(auto i:final) FASTIO::writeInt(i,' ');
		FASTIO::writeChar('\n');
	}
}