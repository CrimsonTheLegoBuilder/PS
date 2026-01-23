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
        double N = scan.nextDouble();
        N /= Math.PI;
        //System.out.println(Math.sqrt(N) * 2);
        System.out.printf("%.15f\n", Math.sqrt(N) * 2);
    }
}