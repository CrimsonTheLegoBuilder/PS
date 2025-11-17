#import<bits/stdc++.h>
using namespace std;
typedef pair<long,long> pp;
long m,n,p1,p2,q,v=2;
pp a[60006],b[300006],z[300006];
int C(pp a,pp b){return a>b;}
int ccw(pp a,pp b,pp c)
{
	long z=(b.first-a.first)*(c.second-a.second)-(b.second-a.second)*(c.first-a.first);
	return z>0?1:z<0?-1:0; // 1 : L // 0 : F // -1 : R
}
int main()
{
	ios::sync_with_stdio(0);
	cin.tie(0);
	cin>>n>>m>>q;
	p1=n;
	for(int i=n;i;i--)
    {
        cin>>a[i].first>>a[i].second;
        if(a[i]<a[p1])p1=i;
    }
    for(int i=1;i<p1;i++)a[n+i]=a[i];
    for(int i=1;i<=n;i++)a[i]=a[i+p1-1];
    for(int i=1;i<=n;i++)
    {
        a[n+i]=a[i];
        if(a[i]>a[v])v=i;
    }
    for(int i=1;i<=m;i++)
    {
        cin>>b[i].first>>b[i].second;
        int c1=ccw(a[1],a[v],b[i]),c2;
        if(c1==1)
        {
            int bx=1,by,bz=v-1;
            for(;bx<bz;)
            {
                by=bx+bz+1>>1;
                if(a[by].first<=b[i].first)bx=by;
                else bz=by-1;
            }
            c2=ccw(a[bx],a[bx+1],b[i]);
            if(c2==-1)z[i]={-1,-1};
            else
            {
                p1=bx;
                if(b[i].first<=a[1].first)p2=v;
                else p2=n;
            }
        }
        if(c1==0)
        {
            if(b[i].first<a[1].first)
            {
                p1=1;
                p2=max(2L,v-1);
            }
            else if(b[i].first>a[v].first)
            {
                p1=v;
                if(v==n)p2=1;
                else p2=n;
            }
            else z[i]={-1,-1};
        }
        if(c1==-1)
        {
            int bx=v,by,bz=n;
            for(;bx<bz;)
            {
                by=bx+bz+1>>1;
                if(a[by].first>=b[i].first)bx=by;
                else bz=by-1;
            }
            c2=ccw(a[bx],a[bx+1],b[i]);
            if(c2==-1)z[i]={-1,-1};
            else
            {
                p1=bx;
                if(b[i].first<=a[1].first)p2=v-1;
                else p2=1;
            }
        }
        if(!z[i].first)
        {
            int bx=p2,by,bz=p1;
            if(bx>bz)bz+=n;
            for(;bx<bz;)
            {
                by=bx+bz>>1;
                if(ccw(a[by],a[by+1],b[i])==-1)bx=by+1;
                else bz=by;
            }
            if(bx>n)bx-=n;
            z[i].first=bx;
            bx=p1,bz=p2;
            if(bx>bz)bz+=n;
            for(;bx<bz;)
            {
                by=bx+bz+1>>1;
                if(ccw(a[by],a[by+1],b[i])==-1)bz=by-1;
                else bx=by;
            }
            if(bx>n)bx-=n;
            z[i].second=bx;
        }
    }
    for(;q--;)
    {
        cin>>p1>>p2;
        if(z[p1].first==-1&&z[p2].first==-1)cout<<n<<"\n";
        else if(z[p1].first==-1)
        {
            int k=z[p2].first-z[p2].second;
            if(k>0)k-=n;
            cout<<n+k+1<<"\n";
        }
        else if(z[p2].first==-1)
        {
            int k=z[p1].first-z[p1].second;
            if(k>0)k-=n;
            cout<<n+k+1<<"\n";
        }
        else
        {
            int c1=z[p1].first;
            int c2=z[p1].second-c1+2,c3=z[p2].first-c1+1,c4=z[p2].second-c1+2;
            if(c2<1)c2+=n;
            if(c2>n)c2-=n;
            if(c3<1)c3+=n;
            if(c4<1)c4+=n;
            if(c4>n)c4-=n;
            if(c3<c4)
            {
                if(c2<c3)cout<<n-c2+c3-c4+5<<"\n";
                else if(c2==c3)
                {
                    if(ccw(b[p1],a[z[p2].first],b[p2])==-1)cout<<n-c2+c3-c4+5<<"\n";
                    else cout<<n-c2+c3-c4+4<<"\n";
                }
                else
                {
                    c3=z[c2>c4?p1:p2].second+1;
                    if(ccw(a[c1],b[p1],b[p2])==-1&&ccw(b[p1],b[p2],a[c3])==-1||ccw(a[c1],b[p2],b[p1])==-1&&ccw(b[p2],b[p1],a[c3])==-1)cout<<n-max(c2,c4)+4<<"\n";
                    else cout<<n-max(c2,c4)+3<<"\n";
                }
            }
            else
            {
                int c1=z[p2].first;
                int c2=z[p2].second-c1+2,c3=z[p1].first-c1+1,c4=z[p1].second-c1+2;
                if(c2<1)c2+=n;
                if(c2>n)c2-=n;
                if(c3<1)c3+=n;
                if(c4<1)c4+=n;
                if(c4>n)c4-=n;
                if(c3<c4)
                {
                    if(c2==c3)
                    {
                        if(ccw(b[p2],a[z[p1].first],b[p1])==-1)cout<<n-c2+c3-c4+5<<"\n";
                        else cout<<n-c2+c3-c4+4<<"\n";
                    }
                    else
                    {
                        c3=z[c2>c4?p2:p1].second+1;
                        if(ccw(a[c1],b[p1],b[p2])==-1&&ccw(b[p1],b[p2],a[c3])==-1||ccw(a[c1],b[p2],b[p1])==-1&&ccw(b[p2],b[p1],a[c3])==-1)cout<<n-max(c2,c4)+4<<"\n";
                        else cout<<n-max(c2,c4)+3<<"\n";
                    }
                }
                else if(1==c4&&c2==c3)
                {
                    if(ccw(b[p1],a[z[p1].first],b[p2])&&ccw(b[p1],a[z[p2].first],b[p2]))cout<<"4\n";
                    else cout<<"3\n";
                }
                else cout<<"3\n";
            }
        }
    }
}