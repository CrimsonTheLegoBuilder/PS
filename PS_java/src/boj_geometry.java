import java.util.Scanner;
import java.util.ArrayList;
import java.util.Collections;

class Pos implements Comparable<Pos> {
    long x, y;
    public static ArrayList<Pos> P = new ArrayList<>();
    public Pos(long x, long y) { this.x = x; this.y = y; }
    @Override
    public int compareTo(Pos rhs) {
        if (this.x != rhs.x) return Long.compare(this.x, rhs.x);
        return Long.compare(this.y, rhs.y);
    }
    public static Pos add(Pos p1, Pos p2) {
        return new Pos(p1.x + p2.x, p1.y + p2.y);
    }
    public static Pos sub(Pos p1, Pos p2) { return new Pos(p1.x - p2.x, p1.y - p2.y); }
    public static long cross(Pos p1, Pos p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
    public static long dot(Pos p1, Pos p2) { return p1.x * p2.x + p1.y * p2.y; }
    public Pos mul(long n) { return new Pos(x * n, y * n); }
    public Pos div(long n) { return new Pos(x / n, y / n); }
    public long Euc() { return x * x + y * y; }
    public double mag() { return Math.sqrt(Euc()); }
    public static long cross3(Pos p1, Pos p2, Pos p3) { return cross(sub(p2, p1), sub(p3, p2)); }
    public static long dot3(Pos p1, Pos p2, Pos p3) { return dot(sub(p2, p1), sub(p3, p2)); }
//    public static Pos intersection(Pos p1, Pos p2, Pos q1, Pos q2) {
//        double a1 = cross3(q1, q2, p1);
//        double a2 = -cross3(q1, q2, p2);
//        return add(p1.mul(a2), p2.mul(a1)).div(a1 + a2);
//    }
}
//public class Main {
public class boj_geometry {
    public static void DEBUG() { System.out.println("Geometry Fuck!!"); }
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        long x = scan.nextLong();
        long y = scan.nextLong();
        long x1 = scan.nextLong();
        long y1 = scan.nextLong();
        long x2 = scan.nextLong();
        long y2 = scan.nextLong();
        Pos p0 = new Pos(x, y);
        Pos p1 = new Pos(x1, y1);
        Pos p2 = new Pos(x2, y1);
        Pos p3 = new Pos(x2, y2);
        Pos p4 = new Pos(x1, y2);

        ArrayList<Double> VD = new ArrayList<>();

        double d1 = Math.abs(Pos.cross3(p1, p2, p0) / Pos.sub(p1, p2).mag());
        if (Pos.dot3(p1, p2, p0) <= 0 && Pos.dot3(p2, p1, p0) <= 0) VD.add(d1);
        double d2 = Math.abs(Pos.cross3(p2, p3, p0) / Pos.sub(p2, p3).mag());
        if (Pos.dot3(p2, p3, p0) <= 0 && Pos.dot3(p3, p2, p0) <= 0) VD.add(d2);
        double d3 = Math.abs(Pos.cross3(p3, p4, p0) / Pos.sub(p3, p4).mag());
        if (Pos.dot3(p3, p4, p0) <= 0 && Pos.dot3(p4, p3, p0) <= 0) VD.add(d3);
        double d4 = Math.abs(Pos.cross3(p4, p1, p0) / Pos.sub(p4, p1).mag());
        if (Pos.dot3(p4, p1, p0) <= 0 && Pos.dot3(p1, p4, p0) <= 0) VD.add(d4);
        double d5 = Pos.sub(p0, p1).mag(); VD.add(d5);
        double d6 = Pos.sub(p0, p2).mag(); VD.add(d6);
        double d7 = Pos.sub(p0, p3).mag(); VD.add(d7);
        double d8 = Pos.sub(p0, p4).mag(); VD.add(d8);
        Collections.sort(VD);
        System.out.printf("%.4f\n", VD.getFirst());
    }
}