import java.io.*;
import java.util.*;

public class boj_bronze {
//public class Main {
    static int N, M;
    static char[][] G;
    static boolean[][] V;
    static int[] dr = { -1, 1, 0, 0, 1, -1, -1, 1 };
    static int[] dc = { 0, 0, -1, 1, 1, 1, -1, -1 };
    public static void main(String[] args) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        //StringBuilder sb = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
        StringTokenizer st = new StringTokenizer(br.readLine());
        N = Integer.parseInt(st.nextToken());
        for (int i = 1; i <= N; i++) {
            String l = br.readLine();
            String[] W = l.split(" ");
            Stack<String> S = new Stack<>();
            for (String w : W) {
                S.push(w);
            }
            bw.write("Case #" + i + ": ");
            while (!S.isEmpty()) {
                bw.write(S.pop());
                if (!S.isEmpty()) bw.write(" ");
            }
            bw.newLine();
        }

        //String ret = String.valueOf(sp + 1);
        //bw.write(ret);
        //bw.newLine();
        bw.flush();
        bw.close();
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
}


/*

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