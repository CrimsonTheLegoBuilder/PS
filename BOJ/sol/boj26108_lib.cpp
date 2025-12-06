#include <bits/stdc++.h>
#define fw fwrite(obuf,p3-obuf,1,stdout)
#define getchar() (p1==p2&&(p2=(p1=buf)+fread(buf,1,1<<20,stdin),p1==p2)?EOF:*p1++)
#define putchar(x) (p3-obuf<1<<20?(*p3++=(x)):(fw,p3=obuf,*p3++=(x)))
#define il inline
using namespace std;

char buf[1 << 20], obuf[1 << 20], *p1 = buf, *p2 = buf, *p3 = obuf, str[20 << 2];
il int read() {
    int x = 0, f = 1;
    char ch = getchar();

    while (!isdigit(ch))
        ch == '-' && (f = -1), ch = getchar();

    while (isdigit(ch))
        x = (x << 3) + (x << 1) + (ch ^ 48), ch = getchar();

    return x * f;
}
template<typename T>
il void write(T x, char sf = '\n') {
    if (x < 0)
        putchar('-'), x = ~x + 1;

    int top = 0;

    do
        str[top++] = x % 10, x /= 10;

    while (x);

    while (top)
        putchar(str[--top] + 48);

    if (sf^'#')
        putchar(sf);
}
using ll = long long;
constexpr int MAXN = 1e5 + 5, MAXK = 105;
int n, m, id[MAXN], pos[MAXN], stk[MAXN];
struct P {
    int x, y;
    il bool operator<(const P &b)const {
        return x != b.x ? x < b.x : y < b.y;
    }
    friend il P operator-(P a, P b) {
        return {a.x - b.x, a.y - b.y};
    }
    friend il ll operator*(P a, P b) {
        return 1ll * a.x * b.y - 1ll * a.y * b.x;
    }
} tp1[MAXN], tp2[MAXN];
struct {
#define lp st[p].lc
#define rp st[p].rc
#define lq st[q].lc
#define rq st[q].rc
    int cur, tot, rt[MAXK];
    int col[MAXN];
    struct SegTree {
        int lc, rc;
        int mn, mx;
        ll ans;
    } st[MAXN << 5];
    P sht[MAXN];
    vector<int>v[MAXK];
    il int build(int x, int s, int t) {
        int p = ++tot;
        st[p] = {0, 0, x, x, 0};

        if (s == t)
            return p;

        int mid = (s + t) >> 1;

        if (x <= mid)
            lp = build(x, s, mid);
        else
            rp = build(x, mid + 1, t);

        return p;
    }
    il int psu(int p) {
        if (!lp || !rp) {
            int t = lp | rp;
            st[p].mn = st[t].mn, st[p].mx = st[t].mx, st[p].ans = st[t].ans;
        } else {
            st[p].mn = st[lp].mn;
            st[p].mx = st[rp].mx;
            P a = sht[st[lp].mn], b = sht[st[lp].mx], c = sht[st[rp].mn], d = sht[st[rp].mx];
            st[p].ans = st[lp].ans + st[rp].ans + (c - a) * (b - a) + (d - a) * (c - a);
        }

        return p;
    }
    il int mge(int l, int r, int p, int q) {
        if (!p || !q)
            return p | q;

        int id = ++tot, mid = (l + r) >> 1;
        st[id].lc = mge(l, mid, lp, lq);
        st[id].rc = mge(mid + 1, r, rp, rq);
        return psu(id);
    }
    il int split(int l, int r, int s, int t, int p) {
        if (!p || l > st[p].mx || r < st[p].mn || l > r)
            return 0;

        if (l <= s && t <= r)
            return p;

        int mid = (s + t) >> 1;
        int x = split(l, r, s, mid, lp), y = split(l, r, mid + 1, t, rp);

        if (!x && !y)
            return 0;

        int res = ++tot;
        st[res].lc = x, st[res].rc = y;
        return psu(res);
    }
    il int askl(int x, int s, int t, int p) {
        if (!p || x < st[p].mn)
            return 0;

        if (s == t)
            return s;

        int mid = (s + t) >> 1;
        P a = sht[st[lp].mx], b = sht[st[rp].mn], c = sht[x];

        if (!rp || x < st[rp].mn || (lp && (c - b) * (b - a) <= 0))
            return askl(x, s, mid, lp);
        else
            return askl(x, mid + 1, t, rp);
    }
    il int askr(int x, int s, int t, int p) {
        if (!p || x > st[p].mx)
            return n + 1;

        if (s == t)
            return s;

        int mid = (s + t) >> 1;
        P a = sht[x], b = sht[st[lp].mx], c = sht[st[rp].mn];

        if (!lp || x > st[lp].mx || (rp && (c - b) * (b - a) <= 0))
            return askr(x, mid + 1, t, rp);
        else
            return askr(x, s, mid, lp);
    }
    il int ins(int A, int B, int l, int r) {
        B = split(l, r, 1, n, B);

        if (!A || !B)
            return A | B;

        int L = askl(st[B].mn, 1, n, A), R = askr(st[B].mx, 1, n, A);
        int C = split(1, L, 1, n, A), D = split(R, n, 1, n, A);
        return mge(1, n, mge(1, n, C, B), D);
    }
    il void init() {
        memset(col, -1, (n + 1) << 2);

        for (int d = 0; d <= 100; d++) {
            int top = 0;

            for (int i = 1; i <= n; i++) {
                if (~col[i])
                    continue;

                while (top > 1 && (sht[i] - sht[stk[top]]) * (sht[stk[top]] - sht[stk[top - 1]]) <= 0)
                    top--;

                stk[++top] = i;
            }

            for (int i = 1; i <= top; i++) {
                col[stk[i]] = d;
                rt[d] = mge(1, n, rt[d], build(stk[i], 1, n));
            }
        }

        cur = tot;
    }
    il ll query(const vector<int> &t) {
        tot = cur;

        for (int i = 0; i <= 100; i++)
            v[i].clear();

        for (auto p : t)
            v[col[p]].push_back(p);

        int i = 0;

        while (!v[i].empty())
            i++;

        int res = rt[i];

        while (i--) {
            int l = 0;
            sort(v[i].begin(), v[i].end());

            for (auto p : v[i])
                res = ins(res, rt[i], l + 1, p - 1), l = p;

            res = ins(res, rt[i], l + 1, n);
        }

        return st[res].ans;
    }
} T1, T2;

int main() {
    n = read(), m = read();

    for (int i = 1; i <= n; i++) {
        int x = read(), y = read();
        T1.sht[i] = {x, y}, T2.sht[i] = {x, -y};
        id[i] = i;
    }

    sort(id + 1, id + n + 1, [](int x, int y) {
        return T1.sht[x] < T1.sht[y];
    });

    for (int i = 1; i <= n; i++)
        pos[id[i]] = i;

    for (int i = 1; i <= n; i++)
        tp1[i] = T1.sht[id[i]], tp2[i] = T2.sht[id[i]];

    swap(T1.sht, tp1), swap(T2.sht, tp2);
    T1.init(), T2.init();
    ll lst = n - 1;

    while (m--) {
        int k = read();
        vector<int>d;

        for (int i = 1; i <= k; i++) {
            ll c = read();
            c = (c + lst) % n + 1;
            d.push_back(pos[c]);
        }

        write(lst = T1.query(d) + T2.query(d));
    }

    return fw, 0;
}