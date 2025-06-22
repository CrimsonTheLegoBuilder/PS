#include<bits/stdc++.h>
using namespace std;

const int N=100100;
const int M=1e6;
typedef long long ll;

struct Point{
    ll x,y,z,r;
    Point(ll x=0,ll y=0,ll z=0,ll r=0):x(x),y(y),z(z),r(r){}
    Point operator * (const int&a)const{return Point(x*a,y*a,z*a,r*a);}
    Point operator - (const Point&a)const{return Point(x-a.x,y-a.y,z-a.z,r+a.r);}

    bool chk_dis(){
        return x*x+y*y+z*z-r*r<=0;
    }
    int get_dis(){
        return max((ll)ceil((sqrt(x*x+y*y+z*z)-r)/2),0ll);
    }
    bool operator < (const Point&a)const{
        return r<a.r;
    }
}a[N];

vector<pair<int,int> >S;

unordered_map<ll,vector<int> >mp(1<<20);

int K;
ll hsh(int x,int y,int z){
    if(x<0||x>M||y<0||y>M||z<0||z>M)return -1;
    return (ll)x*(M+1)*(M+1)+(ll)y*(M+1)+z;
}



void solve1(int n){
    if(!n)
        return;
    int lim=(a[n].r+1)/2;
    int m=1;
    while(a[m].r<lim)
        ++m;
    mp.clear();
    lim*=4;
    for(int i=n;i>=m;--i){
        int x=a[i].x/lim,y=a[i].y/lim,z=a[i].z/lim;
        for(int ix=-1;ix<=1;++ix){
            for(int iy=-1;iy<=1;++iy){
                for(int iz=-1;iz<=1;++iz){
                    auto it=mp.find(hsh(x+ix,y+iy,z+iz));
                    if(it==mp.end())
                        continue;
                    for(auto t:it->second)
                        if((a[i]-a[t]).chk_dis()){
                            S.emplace_back(i,t);
                            if(S.size()>=K)
                                return;
                        }
                }
            }
        }
        mp[hsh(x,y,z)].push_back(i);
    }
    for(int i=m-1;i>=1;--i){
        int x=a[i].x/lim,y=a[i].y/lim,z=a[i].z/lim;
        for(int ix=-1;ix<=1;++ix){
            for(int iy=-1;iy<=1;++iy){
                for(int iz=-1;iz<=1;++iz){
                    auto it=mp.find(hsh(x+ix,y+iy,z+iz));
                    if(it==mp.end())
                        continue;
                    for(auto t:it->second)
                        if((a[i]-a[t]).chk_dis()){
                            S.emplace_back(i,t);
                            if(S.size()>=K)
                                return;
                        }
                }
            }
        }
    }
    solve1(m-1);
}

int n;

bool chk(int m){
    for(int i=1;i<=n;++i)
        a[i].r+=m;
    S.clear();
    solve1(n);
    for(int i=1;i<=n;++i)
        a[i].r-=m;
    return S.size()>=K;
}

vector<ll>ans;
void work(int m){
    S.clear();
    if(m){
        for(int i=1;i<=n;++i)
            a[i].r+=m-1;
        solve1(n);
        for(int i=1;i<=n;++i)
            a[i].r-=m-1;
    }
    ans.clear();
    for(auto i:S)
        ans.push_back((a[i.first]-a[i.second]).get_dis());
    while(ans.size()<K)
        ans.push_back(m);
    sort(ans.begin(),ans.end());
    for(auto i:ans)
        cout<<i<<'\n';
}

int T;
int main(){
    ios::sync_with_stdio(0);cin.tie(0);
    for(cin>>T;T --> 0;){
        cin>>n>>K;
        for(int i=1;i<=n;++i)
            cin>>a[i].x>>a[i].y>>a[i].z>>a[i].r,a[i]=a[i]*2;
        sort(a+1,a+n+1);
        int l=0,r=2e6;
        while(l<r){
            int mid=(l+r)/2;
            if(chk(mid))
                r=mid;
            else
                l=mid+1;
        }
        work(l);
    }
    return 0;
}