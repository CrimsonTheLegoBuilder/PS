#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
typedef long long ll;
typedef long double ld;
const int LEN = 205;
ll sq(const ll& x) { return x * x; }

ll R_, L, X[LEN], Y[LEN], R[LEN];
int N, M, W[LEN];
bool V[LEN];
std::vector<int> G[LEN];
std::queue<int> Q;
int bfs(int c_ = 0) {
    int c = c_;
    Q.push(0);
    while (Q.size()) {
        int p = Q.front(); Q.pop();
        if (V[p]) continue;
        V[p] = 1;
        c += W[p];
        for (int w : G[p]) {
            if (!V[w]) {
                Q.push(w);
            }
        }
    }
    return c;
}
void solve() {
    std::cin.tie(0)->sync_with_stdio(0);
    std::cout.tie(0);
    //freopen("../tests/duck/in_re/82.in", "r", stdin);
    //freopen("../tests/duck/out/82.out", "w", stdout);
    std::cin >> R_ >> L;
    std::cin >> N >> M;
    for (int i = 1; i <= N; i++) std::cin >> X[i] >> Y[i] >> R[i];
    for (int i = N + 1; i <= N + M; i++) std::cin >> X[i] >> Y[i], W[i] = 1;
    for (int i = 1; i <= N; i++) {
        ll D = sq(R_ - R[i]);
        ll d = sq(X[i]) + sq(Y[i]);
        if (R_ <= R[i] || D <= d) {
            G[0].push_back(i);
            G[i].push_back(0);
        }
        for (int j = i + 1; j <= N; j++) {
            D = sq(R[i] + R[j]);
            d = sq(X[i] - X[j]) + sq(Y[i] - Y[j]);
            if (D >= d) {
                G[i].push_back(j);
                G[j].push_back(i);
            }
        }
    }
    ll save = sq(std::max(R_ - L, 0ll));
    int c_ = 0;
    for (int i = N + 1; i <= N + M; i++) {
        ll d = sq(X[i]) + sq(Y[i]);
        if (save <= d) {
            c_ += W[i];
            continue;
        }
        for (int j = 1; j <= N; j++) {
            d = sq(X[i] - X[j]) + sq(Y[i] - Y[j]);
            ll D = sq(R[j] + L);
            if (D >= d) {
                G[j].push_back(i);
            }
        }
    }
    //std::cout << c_ << "\n";
    int total = bfs(c_);
    std::cout << total << "\n";
    return;
}
int main() { solve(); return 0; }
