#include <bits/stdc++.h>
using namespace std;
using ll = long long;

struct Rect {
    ll x1, y1, x2, y2;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int W, H, N;
    cin >> W >> H >> N;
    vector<Rect> ads(N);
    for (int i = 0; i < N; i++) {
        cin >> ads[i].x1 >> ads[i].y1 >> ads[i].x2 >> ads[i].y2;
    }

    // --- 1. 좌표 압축 ---
    vector<ll> xs = { 0, W };
    vector<ll> ys = { 0, H };
    for (auto& r : ads) {
        xs.push_back(r.x1);
        xs.push_back(r.x2);
        ys.push_back(r.y1);
        ys.push_back(r.y2);
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    int XN = xs.size();
    int YN = ys.size();

    // --- 2. 격자 생성 ---
    vector<vector<bool>> blocked(YN - 1, vector<bool>(XN - 1, false));

    for (auto& r : ads) {
        int x1 = lower_bound(xs.begin(), xs.end(), r.x1) - xs.begin();
        int x2 = lower_bound(xs.begin(), xs.end(), r.x2) - xs.begin();
        int y1 = lower_bound(ys.begin(), ys.end(), r.y1) - ys.begin();
        int y2 = lower_bound(ys.begin(), ys.end(), r.y2) - ys.begin();
        for (int y = y1; y < y2; y++) {
            for (int x = x1; x < x2; x++) {
                blocked[y][x] = true;
            }
        }
    }

    // --- 3. 최대 직사각형 찾기 ---
    ll ans = 0;
    vector<ll> height(XN - 1, 0);

    for (int yb = 0; yb < YN - 1; yb++) {
        fill(height.begin(), height.end(), 0);
        for (int yt = yb; yt < YN - 1; yt++) {
            ll rowHeight = ys[yt + 1] - ys[yt];
            for (int x = 0; x < XN - 1; x++) {
                if (blocked[yt][x]) height[x] = 0;
                else height[x] += rowHeight;
            }

            // Largest Rectangle in Histogram (실제 폭 적용)
            vector<int> st;
            vector<ll> left(XN - 1), right(XN - 1);
            for (int i = 0; i < XN - 1; i++) {
                while (!st.empty() && height[st.back()] >= height[i]) st.pop_back();
                left[i] = st.empty() ? 0 : st.back() + 1;
                st.push_back(i);
            }
            st.clear();
            for (int i = XN - 2; i >= 0; i--) {
                while (!st.empty() && height[st.back()] >= height[i]) st.pop_back();
                right[i] = st.empty() ? XN - 2 : st.back() - 1;
                st.push_back(i);
            }
            for (int i = 0; i < XN - 1; i++) {
                if (height[i] == 0) continue;
                ll width = xs[right[i] + 1] - xs[left[i]];
                ans = max(ans, height[i] * width);
            }
        }
    }

    cout << ans << "\n";
    return 0;
}
