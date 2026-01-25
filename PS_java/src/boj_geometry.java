import java.util.Scanner;
import java.util.ArrayList;
import java.util.Collections;

class Pos implements Comparable<Pos> {
    double x, y;
    public static ArrayList<Pos> P = new ArrayList<>();
    public Pos(double x, double y) { this.x = x; this.y = y; }
    @Override
    public int compareTo(Pos rhs) {
        if (this.x != rhs.x) return Double.compare(this.x, rhs.x);
        return Double.compare(this.y, rhs.y);
    }
    public static Pos add(Pos p1, Pos p2) {
        return new Pos(p1.x + p2.x, p1.y + p2.y);
    }
    public static Pos sub(Pos p1, Pos p2) { return new Pos(p1.x - p2.x, p1.y - p2.y); }
    public static double cross(Pos p1, Pos p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
    public static double dot(Pos p1, Pos p2) { return p1.x * p2.x + p1.y * p2.y; }
    public Pos mul(double n) { return new Pos(x * n, y * n); }
    public Pos div(double n) { return new Pos(x / n, y / n); }
    public static double cross3(Pos p1, Pos p2, Pos p3) { return cross(sub(p2, p1), sub(p3, p2)); }
    public static double dot(Pos p1, Pos p2, Pos p3) { return dot(sub(p2, p1), sub(p3, p2)); }
    public static Pos intersection(Pos p1, Pos p2, Pos q1, Pos q2) {
        double a1 = cross3(q1, q2, p1);
        double a2 = -cross3(q1, q2, p2);
        return add(p1.mul(a2), p2.mul(a1)).div(a1 + a2);
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