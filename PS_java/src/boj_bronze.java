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
        long a = sc.nextLong();
        long b = sc.nextLong();
        long c = sc.nextLong();
        for (long x = 1; x <= 100; x++) {
            long l = a * a * a;
            long r = a * (b * b + c * c + x * x) + 2 * b * c * x;
            if (l == r) {
                System.out.println(x);
                return;
            }
        }
        System.out.println(-1);
        //System.out.printf("%.15f\n", Math.max(d1, d2));
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