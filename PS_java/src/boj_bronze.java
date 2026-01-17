//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;

public class boj_bronze {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        int aa = sc.nextInt();
        int ah = sc.nextInt();
        int ba = sc.nextInt();
        int bh = sc.nextInt();
        int ac = ah / ba;
        int bc = bh / aa;
        int c = Math.min(ac, bc);
        ah -= c * ba;
        bh -= c * aa;
        if (ah <= 0 && bh <= 0) {
            System.out.println("DRAW");
            return;
        }
        if (ah <= 0) {
            System.out.println("PLAYER B");
            return;
        }
        if (bh <= 0) {
            System.out.println("PLAYER A");
            return;
        }
        ah -= ba; bh -= aa;
        if (ah <= 0 && bh <= 0) {
            System.out.println("DRAW");
            return;
        }
        if (ah <= 0) {
            System.out.println("PLAYER B");
            return;
        }
        if (bh <= 0) {
            System.out.println("PLAYER A");
            return;
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