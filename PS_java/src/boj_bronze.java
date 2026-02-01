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

        StringTokenizer st = new StringTokenizer(br.readLine());
        int N = Integer.parseInt(st.nextToken());
        N++;
        long[] H = new long[N];
        st = new StringTokenizer(br.readLine());
        for (int i = 0; i < N; i++) {
            H[i] = Long.parseLong(st.nextToken());
        }
        long[] W = new long[N];
        st = new StringTokenizer(br.readLine());
        for (int i = 0; i < N - 1; i++) {
            W[i] = Long.parseLong(st.nextToken());
        }
        double T = 0.;
        long S = 0;
        for (int i = 0; i < N - 1; i++) {
            long h1 = H[i], h2 = H[i + 1];
            long w = W[i];
            S += (h1 + h2) * w;
        }
        T = S * .5;
        String ret = String.format("%.9f", T);
        //String ret = String.valueOf(T);
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