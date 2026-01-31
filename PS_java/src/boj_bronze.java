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
        int T = Integer.parseInt(st.nextToken());
        long s = 0;
        for (int t = 0; t < T; t++) {
            st = new StringTokenizer(br.readLine());
            long x = Long.parseLong(st.nextToken());
            long y = Long.parseLong(st.nextToken());
            long d = x * x + y * y;
            if (d <= (10 * 10)) s += 10;
            else if (d <= (30 * 30)) s += 9;
            else if (d <= (50 * 50)) s += 8;
            else if (d <= (70 * 70)) s += 7;
            else if (d <= (90 * 90)) s += 6;
            else if (d <= (110 * 110)) s += 5;
            else if (d <= (130 * 130)) s += 4;
            else if (d <= (150 * 150)) s += 3;
            else if (d <= (170 * 170)) s += 2;
            else if (d <= (190 * 190)) s += 1;
        }
        //String ret = String.format("%.9f", ans);
        String ret = String.valueOf(s);
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