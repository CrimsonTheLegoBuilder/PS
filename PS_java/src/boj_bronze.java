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
        for (int i = 1; i <= N; i++) {
            int[] A = new int[3];
            for (int j = 0; j < 3; j++) A[j] = sc.nextInt();
            Arrays.sort(A);
            int d = A[0] * A[0] + A[1] * A[1];
            if (d == A[2] * A[2]) System.out.println("Case #" + i + ": YES");
            else System.out.println("Case #" + i + ": NO");
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