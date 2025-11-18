#ifndef LOCAL
#pragma GCC optimize ("Ofast")
#pragma GCC optimize ("unroll-loops")
#endif

#include <bits/stdc++.h>
using namespace std;

using ll=long long;
#define int ll

#define rng(i,a,b) for(int i=int(a);i<int(b);i++)
#define rep(i,b) rng(i,0,b)
#define gnr(i,a,b) for(int i=int(b)-1;i>=int(a);i--)
#define per(i,b) gnr(i,0,b)
#define pb push_back
#define eb emplace_back
#define a first
#define b second
#define bg begin()
#define ed end()
#define all(x) x.bg,x.ed
#define si(x) int(x.size())
#ifdef LOCAL
#define dmp(x) cerr<<__LINE__<<" "<<#x<<" "<<x<<endl
#else
#define dmp(x) void(0)
#endif

template<class t,class u> bool chmax(t&a,u b){if(a<b){a=b;return true;}else return false;}
template<class t,class u> bool chmin(t&a,u b){if(b<a){a=b;return true;}else return false;}

template<class t> using vc=vector<t>;
template<class t> using vvc=vc<vc<t>>;

using pi=pair<int,int>;
using vi=vc<int>;

template<class t,class u>
ostream& operator<<(ostream& os,const pair<t,u>& p){
	return os<<"{"<<p.a<<","<<p.b<<"}";
}

template<class t> ostream& operator<<(ostream& os,const vc<t>& v){
	os<<"{";
	for(auto e:v)os<<e<<",";
	return os<<"}";
}

#define mp make_pair
#define mt make_tuple
#define one(x) memset(x,-1,sizeof(x))
#define zero(x) memset(x,0,sizeof(x))
#ifdef LOCAL
void dmpr(ostream&os){os<<endl;}
template<class T,class... Args>
void dmpr(ostream&os,const T&t,const Args&... args){
	os<<t<<" ";
	dmpr(os,args...);
}
#define dmp2(...) dmpr(cerr,__LINE__,##__VA_ARGS__)
#else
#define dmp2(...) void(0)
#endif

using uint=unsigned;
using ull=unsigned long long;

template<class t,size_t n>
ostream& operator<<(ostream&os,const array<t,n>&a){
	return os<<vc<t>(all(a));
}

template<int i,class T>
void print_tuple(ostream&,const T&){
}

template<int i,class T,class H,class ...Args>
void print_tuple(ostream&os,const T&t){
	if(i)os<<",";
	os<<get<i>(t);
	print_tuple<i+1,T,Args...>(os,t);
}

template<class ...Args>
ostream& operator<<(ostream&os,const tuple<Args...>&t){
	os<<"{";
	print_tuple<0,tuple<Args...>,Args...>(os,t);
	return os<<"}";
}

ll read(){
	ll i;
	cin>>i;
	return i;
}

vi readvi(int n,int off=0){
	vi v(n);
	rep(i,n)v[i]=read()+off;
	return v;
}

pi readpi(int off=0){
	int a,b;cin>>a>>b;
	return pi(a+off,b+off);
}

template<class t>
void print_single(t x,int suc=1){
	cout<<x;
	if(suc==1)
		cout<<"\n";
	if(suc==2)
		cout<<" ";
}

template<class t,class u>
void print_single(const pair<t,u>&p,int suc=1){
	print_single(p.a,2);
	print_single(p.b,suc);
}

template<class T>
void print_single(const vector<T>&v,int suc=1){
	rep(i,v.size())
		print_single(v[i],i==int(v.size())-1?suc:2);
}

template<class T>
void print_offset(const vector<T>&v,ll off,int suc=1){
	rep(i,v.size())
		print_single(v[i]+off,i==int(v.size())-1?suc:2);
}

template<class T,size_t N>
void print_single(const array<T,N>&v,int suc=1){
	rep(i,N)
		print_single(v[i],i==int(N)-1?suc:2);
}

template<class T>
void print(const T&t){
	print_single(t);
}

template<class T,class ...Args>
void print(const T&t,const Args&...args){
	print_single(t,2);
	print(args...);
}

string readString(){
	string s;
	cin>>s;
	return s;
}

template<class T>
T sq(const T& t){
	return t*t;
}

void YES(bool ex=true){
	cout<<"YES\n";
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
void NO(bool ex=true){
	cout<<"NO\n";
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
void Yes(bool ex=true){
	cout<<"Yes\n";
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
void No(bool ex=true){
	cout<<"No\n";
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
//#define CAPITAL
/*
void yes(bool ex=true){
	#ifdef CAPITAL
	cout<<"YES"<<"\n";
	#else
	cout<<"Yes"<<"\n";
	#endif
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
void no(bool ex=true){
	#ifdef CAPITAL
	cout<<"NO"<<"\n";
	#else
	cout<<"No"<<"\n";
	#endif
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}*/
void possible(bool ex=true){
	#ifdef CAPITAL
	cout<<"POSSIBLE"<<"\n";
	#else
	cout<<"Possible"<<"\n";
	#endif
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}
void impossible(bool ex=true){
	#ifdef CAPITAL
	cout<<"IMPOSSIBLE"<<"\n";
	#else
	cout<<"Impossible"<<"\n";
	#endif
	if(ex)exit(0);
	#ifdef LOCAL
	cout.flush();
	#endif
}

constexpr ll ten(int n){
	return n==0?1:ten(n-1)*10;
}

const ll infLL=LLONG_MAX/3;

#ifdef int
const int inf=infLL;
#else
const int inf=INT_MAX/2-100;
#endif

int topbit(signed t){
	return t==0?-1:31-__builtin_clz(t);
}
int topbit(ll t){
	return t==0?-1:63-__builtin_clzll(t);
}
int topbit(ull t){
	return t==0?-1:63-__builtin_clzll(t);
}
int botbit(signed a){
	return a==0?32:__builtin_ctz(a);
}
int botbit(ll a){
	return a==0?64:__builtin_ctzll(a);
}
int botbit(ull a){
	return a==0?64:__builtin_ctzll(a);
}
int popcount(signed t){
	return __builtin_popcount(t);
}
int popcount(ll t){
	return __builtin_popcountll(t);
}
int popcount(ull t){
	return __builtin_popcountll(t);
}
int bitparity(ll t){
	return __builtin_parityll(t);
}
bool ispow2(int i){
	return i&&(i&-i)==i;
}
ll mask(int i){
	return (ll(1)<<i)-1;
}
ull umask(int i){
	return (ull(1)<<i)-1;
}
ll minp2(ll n){
	if(n<=1)return 1;
	else return ll(1)<<(topbit(n-1)+1);
}

bool inc(int a,int b,int c){
	return a<=b&&b<=c;
}

template<class t> void mkuni(vc<t>&v){
	sort(all(v));
	v.erase(unique(all(v)),v.ed);
}

ll rand_int(ll l, ll r) { //[l, r]
	//#ifdef LOCAL
	static mt19937_64 gen;
	/*#else
	static mt19937_64 gen(chrono::steady_clock::now().time_since_epoch().count());
	#endif*/
	return uniform_int_distribution<ll>(l, r)(gen);
}

ll rand_int(ll k){ //[0,k)
	return rand_int(0,k-1);
}

template<class t>
void myshuffle(vc<t>&a){
	rep(i,si(a))swap(a[i],a[rand_int(0,i)]);
}

template<class t>
int lwb(const vc<t>&v,const t&a){
	return lower_bound(all(v),a)-v.bg;
}
template<class t>
bool bis(const vc<t>&v,const t&a){
	return binary_search(all(v),a);
}

vvc<int> readGraph(int n,int m){
	vvc<int> g(n);
	rep(i,m){
		int a,b;
		cin>>a>>b;
		//sc.read(a,b);
		a--;b--;
		g[a].pb(b);
		g[b].pb(a);
	}
	return g;
}

vvc<int> readTree(int n){
	return readGraph(n,n-1);
}

vc<ll> presum(const vi&a){
	vc<ll> s(si(a)+1);
	rep(i,si(a))s[i+1]=s[i]+a[i];
	return s;
}
//BIT で数列を管理するときに使う (CF850C)
template<class t>
vc<t> predif(vc<t> a){
	gnr(i,1,si(a))a[i]-=a[i-1];
	return a;
}
template<class t>
vvc<ll> imos(const vvc<t>&a){
	int n=si(a),m=si(a[0]);
	vvc<ll> b(n+1,vc<ll>(m+1));
	rep(i,n)rep(j,m)
		b[i+1][j+1]=b[i+1][j]+b[i][j+1]-b[i][j]+a[i][j];
	return b;
}

//verify してないや
void transvvc(int&n,int&m){
	swap(n,m);
}
template<class t,class... Args>
void transvvc(int&n,int&m,vvc<t>&a,Args&...args){
	assert(si(a)==n);
	vvc<t> b(m,vi(n));
	rep(i,n){
		assert(si(a[i])==m);
		rep(j,m)b[j][i]=a[i][j];
	}
	a.swap(b);
	transvvc(n,m,args...);
}
//CF854E
void rotvvc(int&n,int&m){
	swap(n,m);
}
template<class t,class... Args>
void rotvvc(int&n,int&m,vvc<t>&a,Args&...args){
	assert(si(a)==n);
	vvc<t> b(m,vi(n));
	rep(i,n){
		assert(si(a[i])==m);
		rep(j,m)b[m-1-j][i]=a[i][j];
	}
	a.swap(b);
	rotvvc(n,m,args...);
}

//ソートして i 番目が idx[i]
//CF850C
template<class t>
vi sortidx(const vc<t>&a){
	int n=si(a);
	vi idx(n);iota(all(idx),0);
	sort(all(idx),[&](int i,int j){return a[i]<a[j];});
	return idx;
}
//vs[i]=a[idx[i]]
//例えば sortidx で得た idx を使えば単にソート列になって返ってくる
//CF850C
template<class t>
vc<t> a_idx(const vc<t>&a,const vi&idx){
	int n=si(a);
	assert(si(idx)==n);
	vc<t> vs(n);
	rep(i,n)vs[i]=a[idx[i]];
	return vs;
}
//CF850C
vi invperm(const vi&p){
	int n=si(p);
	vi q(n);
	rep(i,n)q[p[i]]=i;
	return q;
}

template<class t,class s=t>
s SUM(const vc<t>&a){
	return accumulate(all(a),s(0));
}

template<class t>
t MAX(const vc<t>&a){
	return *max_element(all(a));
}

template<class t>
t MIN(const vc<t>&a){
	return *min_element(all(a));
}

template<class t>
pair<t,int> MINi(const vc<t>&a){
	auto itr=min_element(all(a));
	return mp(*itr,itr-a.bg);
}

template<class t,class u>
pair<t,u> operator+(const pair<t,u>&a,const pair<t,u>&b){
	return mp(a.a+b.a,a.b+b.b);
}

vi vid(int n){
	vi res(n);iota(all(res),0);
	return res;
}

template<class S>
S getrev(S s){
	reverse(all(s));
	return s;
}

pi operator+(pi a,pi b){return pi(a.a+b.a,a.b+b.b);}

template<class t>
t gpp(vc<t>&vs){
	assert(si(vs));
	t res=move(vs.back());
	vs.pop_back();
	return res;
}

//copied from yosupo's library
//ptARTLY VERIFIED

//USACO 2022 January ptlatinum C

//#define GEOF

#ifdef GEOF
using ld=long double;
//using ld=double;
const ld PI=acos(ld(-1));
#else
using ld=ll;
#endif
const ld eps=1e-9;
int sgn(ld a){return a<-eps?-1:(a>eps?1:0);}
int sgn(ld a,ld b){return sgn(a-b);}
/*
using pt=complex<ld>;
#define x real()
#define y imag()
*/
struct pt {
    ld x,y;
    //pt(ld _x = ld(0), ld _y = ld(0)) : x(_x), y(_y) {}
    pt():x(0),y(0){}
    pt(ld xx,ld yy):x(xx),y(yy){}
    pt operator+(const pt& r) const { return {x + r.x, y + r.y}; }
    pt operator-(const pt& r) const { return {x - r.x, y - r.y}; }
    pt operator*(const pt& r) const {
        return {x * r.x - y * r.y, x * r.y + y * r.x};
    }
    pt inv()const{
		ld d=norm();
		return {x/d,-y/d};
	}
    pt operator/(const pt&r)const{
		return operator*(r.inv());
	}

    pt operator*(const ld& r) const { return {x * r, y * r}; }
    pt operator/(const ld& r) const { return {x / r, y / r}; }

    pt& operator+=(const pt& r) { return *this = *this + r; }
    pt& operator-=(const pt& r) { return *this = *this - r; }
    pt& operator*=(const pt& r) { return *this = *this * r; }
    pt& operator/=(const pt& r) { return *this = *this / r; }
    pt& operator*=(const ld& r) { return *this = *this * r; }
    pt& operator/=(const ld& r) { return *this = *this / r; }

    pt operator-() const { return {-x, -y}; }

	static int cmp(const pt&a,const pt&b){
		int v=sgn(a.x,b.x);
		return v?v:sgn(a.y,b.y);
	}
    bool operator<(const pt& r) const {
        return cmp(*this,r)<0;
    }
    bool operator<=(const pt& r) const {
        return cmp(*this,r)<=0;
    }
    bool operator>(const pt& r) const {
        return cmp(*this,r)>0;
    }
    bool operator==(const pt& r) const { return sgn((*this - r).rabs()) == 0; }
    bool operator!=(const pt& r) const { return !(*this == r); }

    ld norm() const { return x * x + y * y; }
    ld rabs() const { return max(std::abs(x), std::abs(y)); }  // robust abs
    pair<ld, ld> to_pair() const { return {x, y}; }
    #ifdef GEOF
    ld abs() const { return sqrt(norm()); }
    ld arg() const { return atan2(y, x); }
    static pt polar(ld le, ld th) { return {le * cos(th), le * sin(th)}; }
	#endif
};
istream& operator>>(istream& is, pt& p){
	return is>>p.x>>p.y;
}
ostream& operator<<(ostream& os, const pt& p) {
    return os << "pt(" << p.x << ", " << p.y << ")";
}
ld norm(const pt&a){
	return a.norm();
}
#ifdef GEOF
ld abs(const pt&a){
	return sqrt(norm(a));
}
//XXII Opencup Gpt of Ural K
pt normalized(const pt&a){
	return a/abs(a);
}
ld arg(const pt&a){return a.arg();}
//normalize to [-PI,PI)
//Contest 2: ptKU Contest 1, ptTZ Summer 2022 Day 4
ld normarg(ld a){
	ld res=fmod(a+PI,2*PI);
	if(res<0)res+=PI;
	else res-=PI;
	return res;
}
//AOJ1183
//arg between ab
//assume given lengths are valid
ld arg(ld a,ld b,ld c){
	return acos(min(max((a*a+b*b-c*c)/(2*a*b),ld(-1)),ld(1)));
}
#endif
ld norm(const pt&a,const pt&b){
	return (a-b).norm();
}
ld dot(const pt&a,const pt&b){return a.x*b.x+a.y*b.y;}
ld crs(const pt& a,const pt& b){return a.x*b.y-a.y*b.x;}
ld crs(const pt& a,const pt& b,const pt& c){return crs(b-a,c-a);}
int ccw(const pt& a,const pt& b){return sgn(crs(a,b));}
int ccw(const pt& a,const pt& b,const pt& c){return ccw(b-a,c-a);}
//(-pi,0](0,pi]
int argtype(const pt&a){
	if(sgn(a.y)==0)return a.x<0?1:0;
	return a.y<0?0:1;
}
int argcmp(const pt&a,const pt&b){
	int at=argtype(a),bt=argtype(b);
	if(at!=bt)return at<bt?-1:1;
	return -ccw(a,b);
};
bool argless(const pt&a,const pt&b){return argcmp(a,b)<0;}
//(-2)[a,-1](0)[b,1](2)
int bet(pt a,pt b,pt c){
	pt d=b-a;
	ld e=dot(d,c-a);
	if(sgn(e)<=0)return sgn(e)-1;
	return sgn(e-norm(d))+1;
}
//AOJ0153
//三角形 abc に d が含まれるか？0-no,1-edge,2-in
int cont(pt a,pt b,pt c,pt d){
	if(ccw(a,b,c)==-1) swap(b,c);
	return min({ccw(a,b,d),ccw(b,c,d),ccw(c,a,d)})+1;
}
//(a,b) を結ぶ直線を考え，x 座標との交点の y 座標を求める
//(分子,分母)を返す
pair<ld,ld> xcut(const pt&a,const pt&b,ld x){
	return mp(a.y*(b.x-x)-b.y*(a.x-x),b.x-a.x);
}
//XXII Opencup Gpt of Ural K
pt rot90(pt a){
	return pt(-a.y,a.x);
}
#ifdef GEOF
ld xcutval(const pt&a,const pt&b,ld x){
	auto [p,q]=xcut(a,b,x);
	return p/q;
}
//AOJ1183
//Xmas2010 E
//-+ の 順で返す
//a の符号により，small/large が決まる
int qeq(ld a,ld b,ld c,ld&d,ld&e){
	if(sgn(a)==0){
		if(sgn(b)==0)return 0;
		d=-c/b;
		return 1;
	}
	ld f=b*b-4*a*c;
	if(sgn(f)<0)return 0;
	ld g=sqrt(max(f,ld(0)));
	d=(-b-g)/(2*a);
	e=(-b+g)/(2*a);
	return sgn(f)+1;
}
#else
pt normdir(pt a){
	if(a==pt(0,0))return a;
	int g=gcd(a.x,a.y);
	return pt(a.x/g,a.y/g);
}
#endif

ld area2(const vc<pt>&ps){
	ld res=0;
	rep(i,si(ps))res+=crs(ps[i],ps[(i+1)%si(ps)]);
	return res;
}

//The 1st Universal Cup. Stage 15: Hangzhou J
struct S{
	//多分線分同士が端点以外で交わらないというケースだけ動くコード
	//線分に接するみたいな状況が発生すると話が変わってきてしまう
	pt a,b;
	S(){}
	S(const pt&aa,const pt&bb):a(min(aa,bb)),b(max(aa,bb)){assert(a<b);}
	bool operator<(const S&r)const{
		if(a==r.a)return ccw(a,b,r.b)>0;
		else if(a<r.a)return ccw(a,b,r.a)>0;
		else return ccw(r.a,r.b,a)<0;
	}
};

void flip(vc<pt>&vs){
	for(auto&[x,y]:vs)swap(x,y);
	reverse(all(vs));
}

template<class t>
struct BIT{
	vc<t> buf;
	int s;
	BIT(int n=0){init(n);}
	BIT(const vc<t>&a){init(a);}
	void init(int n){buf.clear();buf.resize(s=n);}
	void init(const vc<t>&a){
		s=si(a);
		buf.resize(s);
		rep(i,s)buf[i]=a[i];
		rep(i,s){
			int j=i+((i+1)&(-i-1));
			if(j<s)buf[j]+=buf[i];
		}
	}
	void add(int i,t v){
		for(;i<s;i+=(i+1)&(-i-1))
			buf[i]+=v;
	}
	t get(int i){
		t res=t();
		for(;i>=0;i-=(i+1)&(-i-1))
			res+=buf[i];
		return res;
	}
	t sum(int b,int e){
		return get(e-1)-get(b-1);
	}
	void add_range(int b,int e,t v){
		add(b,v);
		add(e,-v);
	}
	int kth(int k){
		int res=0;
		for(int i=topbit(s);i>=0;i--){
			int w=res+(1<<i);
			if(w<=s&&buf[w-1]<=k){
				k-=buf[w-1];
				res=w;
			}
		}
		return res;
	}
	//yukicoder No.1024
	int kth_helper(int k,int i){
		return kth(k+get(i-1));
	}
};

void slv(){
	int n;cin>>n;
	vc<pt> ps(n);
	rep(i,n)cin>>ps[i];
	
	vc<pt> ans;
	{
		vc<pt> dir;
		rep(i,n)dir.pb(ps[(i+1)%n]-ps[i]);
		sort(all(dir),argless);
		int s=0;
		for(int i=0;i<n;){
			int j=i;
			pt v;
			while(j<n&&ccw(dir[i],dir[j])==0)v+=dir[j++];
			dir[s++]=v;
			i=j;
		}
		dir.resize(s);
		pt cur;
		for(auto v:dir){
			ans.pb(cur);
			cur+=v;
		}
	}

	//여기까진 단순한 각도순 정렬
	
	rep(_,2){
		vc<S> ss(n);
		vc<tuple<pt,int,int>> qs;
		vi kind(n);
		rep(i,n){
			pt a=ps[i],b=ps[(i+1)%n];
			ss[i]=S(a,b);
			qs.eb(ss[i].a,i,0);
			qs.eb(ss[i].b,i,1);
			kind[i]=a<b?0:1;
		}
		sort(all(qs));
		set<pair<S,int>> cur;
		vc<decltype(cur)::iterator> itrs(n);
		BIT<int> bit(n);
		auto work=[&](decltype(cur)::iterator itr,int w){
			if(itr==cur.ed||itr==cur.bg)return;
			int i=prev(itr)->b,j=itr->b;
			if(kind[i]==1&&kind[j]==0){
				if(i<j){
					bit.add_range(0,i+1,w);
					bit.add_range(j+1,n,w);
				}else{
					bit.add_range(j+1,i+1,w);
				}
			}
		};
		int mx=-inf;
		for(auto[x,y]:ps)chmax(mx,x);
		int pre=0;
		for(auto [pos,i,k]:qs){
			if(si(cur)){
				int b=cur.bg->b;
				mx+=(pos.x-pre)*bit.get(b);
			}
			pre=pos.x;
			if(k==0){
				auto obj=mp(ss[i],i);
				auto itr=cur.lower_bound(obj);
				work(itr,-1);
				itr=cur.insert(itr,obj);
				work(itr,1);
				work(next(itr),1);
				itrs[i]=itr;
			}else{
				auto itr=itrs[i];
				work(itr,-1);
				work(next(itr),-1);
				work(cur.erase(itr),1);
			}
		}
		int now=-inf;
		for(auto[x,y]:ans)chmax(now,x);
		mx-=now;
		for(auto&[x,y]:ans)x+=mx;
		flip(ps);
		flip(ans);
	}
	rotate(ans.bg,min_element(all(ans)),ans.ed);
	print(si(ans));
	for(auto[x,y]:ans)print(x,y);
}

signed main(){
	cin.tie(0);
	ios::sync_with_stdio(0);
	cout<<fixed<<setprecision(20);
	
	//int t;cin>>t;rep(_,t)
	slv();
}
