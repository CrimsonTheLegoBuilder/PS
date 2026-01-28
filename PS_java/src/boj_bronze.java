//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;
import java.util.Arrays;

public class boj_bronze {
//public class Main {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        int T = sc.nextInt();
        for (int t = 1; t <= T; t++) {
            int r = sc.nextInt();
            int a = sc.nextInt();
            int b = sc.nextInt();
            double total = r * r * Math.PI;
            int cur = r;
            long cnt = 0;
            while (cur > 0) {
                if (cnt % 2 != 1) cur *= a;
                else cur /= b;
                //System.out.println("DEBUG[" + t + "]:: cur:: " + cur);
                cnt++;
                total += cur * cur * Math.PI;
            }
            System.out.print("Case #" + t + ": ");
            System.out.printf("%.6f\n", total);
        }
    }
}

/*
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