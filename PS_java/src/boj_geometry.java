import java.util.Scanner;
import java.util.ArrayList;
import java.util.Collections;

class Pos implements Comparable<Pos> {
    int x, y;
    public static ArrayList<Pos> P = new ArrayList<>();
    public Pos(int x, int y) { this.x = x; this.y = y; }
    @Override
    public int compareTo(Pos rhs) {
        if (this.x != rhs.x) return this.y - rhs.x;
        return this.y - rhs.y;
    }
    public static Pos sub(Pos p1, Pos p2) {
        return new Pos(p2.x - p1.x, p2.y - p1.y);
    }
    public static long cross(Pos p1, Pos p2) {
        return (long)p1.x * p2.y - (long)p1.y * p2.x;
    }
    public static long cross3(Pos p1, Pos p2, Pos p3) {
        return cross(sub(p1, p2), sub(p2, p3));
    }
}
//public class Main {
public class boj_geometry {
    public static void DEBUG() { System.out.println("Geometry Fuck!!"); }
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        int N = scan.nextInt();
        int x0 = (int)1e9 + 1;
        int x1 = -(int)1e9 - 1;
        int y0 = (int)1e9 + 1;
        int y1 = -(int)1e9 - 1;
        for (int i = 0; i < N; i++) {
            int x = scan.nextInt();
            int y = scan.nextInt();
            x0 = Math.min(x0, x);
            y0 = Math.min(y0, y);
            x1 = Math.max(x1, x);
            y1 = Math.max(y1, y);
        }
        int d = Math.max((x1 - x0), (y1 - y0));
        System.out.println((long)d * d);
    }
}