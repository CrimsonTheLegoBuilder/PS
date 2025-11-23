#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#define right x
#define left y
typedef long long ll;
const int LEN = 2e5 + 1;
const ll MOD = 1'000'000'007;

int N, M, K, S[LEN], E[LEN];// , A[LEN];
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(8);
	int h0, m0, s0, t0 = 0;
	scanf("%d:%d:%d", &h0, &m0, &s0);
	t0 += h0 * 3600;
	t0 += m0 * 60;
	t0 += s0;
	//printf("t0: %d\n", t0);
	int h1, m1, s1, t1 = 0;
	scanf("%d:%d:%d", &h1, &m1, &s1);
	t1 += h1 * 3600;
	t1 += m1 * 60;
	t1 += s1;
	//printf("t1: %d\n", t1);
	int t, k;
	scanf("%d %d", &t, &k);
	N = t1 - t0;
	if (N < 0) { printf("0\n"); return; }
	//printf("N: %d\n", N);
	int ans = 1;
	if (N < (t * (100 - k)) / 100) ans = 0;
	printf("%d\n", ans);
	return;
}
int main() { solve(); return 0; }//boj34702