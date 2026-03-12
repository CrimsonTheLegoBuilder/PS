import java.io.*;
import java.util.*;
import java.math.BigInteger;

public class boj_bronze {
//public class Main {
    static int[] dr = { -1, 1, 0, 0, 1, -1, -1, 1 };
    static int[] dc = { 0, 0, -1, 1, 1, 1, -1, -1 };
    public static void main(String[] args) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        //StringBuilder sb = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
        //StringTokenizer st = new StringTokenizer(br.readLine());
        //st = new StringTokenizer(br.readLine());
        //long a = Long.parseLong(st.nextToken());
        //String A = br.readLine();
        StringTokenizer st;
        st = new StringTokenizer(br.readLine());
        int T = Integer.parseInt(st.nextToken());
        long F = 0;
        int L = 0;
        int R = 0;
        for (int t = 1; t <= T; t++) {
            st = new StringTokenizer(br.readLine());
            int l = Integer.parseInt(st.nextToken());
            int r = Integer.parseInt(st.nextToken());
            if (t >= 2) {
                if (l != 0 && l == L) F++;
                if (r != 0 && r == R) F++;
            }
            if (l != 0 && r != 0 && l == r) F++;
            L = l;
            R = r;
        }
        bw.write(F + "\n");
        bw.flush();
        bw.close();
    }
    static int gcd(int p, int q) {
        if (q == 0) return p;
        return gcd(q, p % q);
    }
}
/*
        //String ret = String.format("%.10f", d);
        //String ret = String.valueOf(ans);
        //bw.write(ans + "\n");
        //bw.newLine();
*/


/*
    static void process(char c) {
        int c1 = 0, c2 = 0;
        switch (c) {
            case 'A': c1 = 1; c2 = 2; break;
            case 'B': c1 = 1; c2 = 3; break;
            case 'C': c1 = 1; c2 = 4; break;
            case 'D': c1 = 2; c2 = 3; break;
            case 'E': c1 = 2; c2 = 4; break;
            case 'F': c1 = 3; c2 = 4; break;
        }
        if (s == c1) s = c2;
        else if (s == c2) s = c1;
        if (l == c1) l = c2;
        else if (l == c2) l = c1;
    }
    static boolean isVowel(char c) {
        c = Character.toLowerCase(c);
        return c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u';
    }
    static void bfs(int r, int c) {
        Queue<int[]> q = new LinkedList<>();
        q.offer(new int[]{ r, c });
        V[r][c] = true;
        while (!q.isEmpty()) {
            int[] cur = q.poll();
            for (int i = 0; i < 4; i++) {
                int nr = cur[0] + dr[i];
                int nc = cur[1] + dc[i];
                if (0 <= nr && nr < N && 0 <= nc && nc < M) {
                    if (G[nr][nc] == '#' && !V[nr][nc]) {
                        V[nr][nc] = true;
                        q.offer(new int[]{ nr, nc });
                    }
                }
            }
        }
    }
    public static void main(String[] args) throws IOException {
        Scanner sc = new Scanner(System.in);
        int T = sc.nextInt();
        for (int t = 1; t <= T; t++) {
            int r = sc.nextInt();
            int b = sc.nextInt();
            double s = (double)r * Math.sqrt(2);
            double ans = s * s;
            if (b < s) {
                int r2 = r * 2;
                double d = Math.sqrt(r2 * r2 - b * b);
                ans = d * b;
            }
            System.out.printf("%.3f\n", ans);
        }
    }

            String s = scan.next();
            int len = s.length();
            char l = s.charAt(len / 2 - 1);
            char r = s.charAt(len / 2);
            if (l == r) {
                System.out.println("Do-it");
            }
            else {
                //System.out.println("Do-it-Not");
            }
 */

        /*
        for (int i = 0; i < 2; i++) {
            st = new StringTokenizer(br.readLine());
            for (int j = 0; j < 3; j++) {
                l[i][j] = Long.parseLong(st.nextToken());
            }
        }
        Arrays.sort(l[0]);
        Arrays.sort(l[1]);
        long ans = 1;
        if (l[0][2] * l[0][2] != l[0][1] * l[0][1] + l[0][0] * l[0][0]) ans = 0;
        if (l[1][2] * l[1][2] != l[1][1] * l[1][1] + l[1][0] * l[1][0]) ans = 0;
        if (l[0][2] != l[1][2] || l[0][1] != l[1][1] || l[0][0] != l[1][0]) ans = 0;
        */