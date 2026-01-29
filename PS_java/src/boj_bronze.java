import java.io.*;
import java.util.Scanner;
import java.util.Arrays;
import java.util.StringTokenizer;

public class boj_bronze {
//public class Main {
    public static void main(String[] args) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        StringBuilder sb = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));

        StringTokenizer st= new StringTokenizer(br.readLine());
        int N = Integer.parseInt(st.nextToken());
        int M = Integer.parseInt(st.nextToken());

        long[][] P = new long[N][2];
        long[][] Q = new long[M][2];

        for (int i = 0; i < N; i++) {
            st = new StringTokenizer(br.readLine());
            P[i][0] = Long.parseLong(st.nextToken());
            P[i][1] = Long.parseLong(st.nextToken());
        }

        for (int j = 0; j < M; j++) {
            st = new StringTokenizer(br.readLine());
            Q[j][0] = Long.parseLong(st.nextToken());
            Q[j][1] = Long.parseLong(st.nextToken());
        }

        long ans = -1;

        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++) {
                long x = Q[j][0] - P[i][0];
                long y = Q[j][1] - P[i][1];
                ans = Math.max(ans, x * x + y * y);
            }
        }

        //String ret = String.format("%.7f", d);
        String ret = String.valueOf(ans);
        bw.write(ret);
        bw.newLine();
        bw.flush();
        bw.close();
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