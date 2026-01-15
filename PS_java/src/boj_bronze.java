//import java.io.BufferedWriter;
//import java.io.OutputStreamWriter;
//import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;

public class boj_bronze {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        int N = sc.nextInt();
        for (int i = 0; i <= N; i++) {
            if (i == 0) { System.out.println("int a;"); }
            else if (i == 1) { System.out.println("int *ptr = &a;"); }
            else {
                System.out.print("int ");
                for (int j = 0; j < i; j++) System.out.print("*");
                if (i - 1 == 1) System.out.println("ptr" + i + " = &ptr;");
                else System.out.println("ptr" + i + " = &ptr" + (i - 1) + ";");
            }
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