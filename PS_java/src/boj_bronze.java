//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;
import java.util.Arrays;

public class boj_bronze {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        int N = sc.nextInt();
        int c = 0;
        for (int i = 1; i <= 5; i++) {
            int j = N - i;
            if (j >= 0 && j <= i) c++;
        }
        System.out.println(c);
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