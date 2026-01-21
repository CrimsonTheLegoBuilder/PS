//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;
import java.util.Arrays;

public class boj_bronze {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        while (true) {
            int c = sc.nextInt();
            int d = sc.nextInt();
            if (c == 0 && d == 0) break;
            int t1 = c * 30 + d * 40;
            int t2 = c * 35 + d * 30;
            int t3 = c * 40 + d * 20;
            System.out.println(Math.min(t1, Math.min(t2, t3)));
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