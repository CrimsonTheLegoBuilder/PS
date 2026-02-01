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
        for (int i = 1; i <= N; i++) {
            int i2 = i * i;
            if (i2 > N) break;
            int j2 = N - i2;
            double j = Math.sqrt(j2);
            int j0 = (int)j;
            int j1 = j0 + 1;

            String s;/* = String.valueOf(i);
            bw.write(s);
            bw.write(" ");
            s = String.valueOf(j0);
            bw.write(s);
            bw.newLine();
            */
            if (i2 + j0 * j0 == N) {
                int x = 0, y = 0;
                int dx = i, dy = j0;
                for (int k = 0; k < 4; k++) {
                    s = String.valueOf(x);
                    bw.write(s);
                    bw.write(" ");
                    s = String.valueOf(y);
                    bw.write(s);
                    bw.newLine();
                    x += dx;
                    y += dy;
                    int tmp = dx;
                    dx = -dy;
                    dy = tmp;
                }
                bw.flush();
                return;
            }
            else if (i2 + j1 * j1 == N) {
                int x = 0, y = 0;
                int dx = i, dy = j1;
                for (int k = 0; k < 4; k++) {
                    s = String.valueOf(x);
                    bw.write(s);
                    bw.write(" ");
                    s = String.valueOf(y);
                    bw.write(s);
                    bw.newLine();
                    x += dx;
                    y += dy;
                    int tmp = dx;
                    dx = -dy;
                    dy = tmp;
                }
                bw.flush();
                return;
            }
        }
        //String ret = String.format("%.1f", r * 2 * Math.PI + .05);
        //String ret = String.valueOf(T);
        bw.write("Impossible");
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