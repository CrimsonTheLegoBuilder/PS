//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;

public class boj_bronze {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        while (true) {
            int l = sc.nextInt();
            int w = sc.nextInt();
            int a = sc.nextInt();
            if (l == 0 && w == 0 && a == 0) break;
            if (l == 0) System.out.println(a / w + " " + w + " " + a);
            if (w == 0) System.out.println(l + " " + a / l + " " + a);
            if (a == 0) System.out.println(l + " " + w + " " + l * w);
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