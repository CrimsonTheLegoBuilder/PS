import java.io.*;
import java.util.*;

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
    public static void main(String[] args) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        //StringBuilder sb = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
        //StringTokenizer st = new StringTokenizer(br.readLine());
        StringTokenizer st;
        st = new StringTokenizer(br.readLine());
        int T = Integer.parseInt(st.nextToken());
        for (int t = 0; t < T; t++) {
            st = new StringTokenizer(br.readLine());
            int N = Integer.parseInt(st.nextToken());
            Pos[] P = new Pos[N];
            Pos a = new Pos(0, 0);
            Pos b = new Pos(0, 0);
            long D = 10000000000000000L;
            for (int i = 0; i < N; i++) {
                st = new StringTokenizer(br.readLine());
                P[i] = new Pos(Long.parseLong(st.nextToken()), Long.parseLong(st.nextToken()));
            }
            Arrays.sort(P);
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    long l = (Pos.sub(P[i], P[j])).Euc();
                    if (D > l) {
                        D = l;
                        a = P[i];
                        b = P[j];
                    }
                }
            }
            //String ret = String.format("%.10f", d);
            //String ret = String.valueOf(ans);
            bw.write(a.x + " " + a.y + " " + b.x + " " + b.y + "\n");
            //bw.newLine();
        }
        bw.flush();
        bw.close();
    }
}