const fs = require('fs');
const assert = require('assert');

class IO {
  constructor() {
    this.stdin = fs.readFileSync(0);
    this.text = this.stdin.toString();
    this.i = 0;
    this.buffer = '';
    process.on('exit', () => this.flush());
  }

  input() {
    while (this.stdin[this.i] <= 32) this.i++;
    const s = this.i;
    while (this.stdin[this.i] > 32) this.i++;
    return this.text.slice(s, this.i);
  }

  eof() {
    while (this.i < this.stdin.length && this.stdin[this.i] <= 32) this.i++;
    return this.i === this.stdin.length;
  }

  inputLine() {
    if (this.stdin[this.i] === 10) this.i++;
    const s = this.i;
    while (this.stdin[this.i] !== 10) this.i++;
    return this.text.slice(s, this.i);
  }

  eofLine() {
    if (this.stdin[this.i] === 10) this.i++;
    return this.i === this.stdin.length;
  }

  inputAll() {
    const s = this.i;
    this.i = this.stdin.length;
    return this.text.slice(s, this.i);
  }

  print(v, end = '\n') {
    this.buffer += v.toString() + end;
    if (this.buffer.length > 1 << 13) this.flush();
  }

  flush() {
    fs.writeSync(1, this.buffer);
    this.buffer = '';
  }
}

const io = new IO();

const eps = 1e-8;

class Point {
  constructor(x = 0, y = 0) {
    this.x = x;
    this.y = y;
  }

  add(p) {
    return new Point(this.x + p.x, this.y + p.y);
  }

  sub(p) {
    return new Point(this.x - p.x, this.y - p.y);
  }

  scalaMul(k) {
    return new Point(this.x * k, this.y * k);
  }

  scalaDiv(k) {
    return new Point(this.x / k, this.y / k);
  }

  gt(p) {
    return this.x !== p.x ? this.x > p.x : this.y > p.y;
  }

  lt(p) {
    return this.x !== p.x ? this.x < p.x : this.y < p.y;
  }

  eq(p) {
    return this.x === p.x && this.y === p.y;
  }

  ge(p) {
    return Math.abs(this.x - p.x) > eps ? this.x + eps >= p.x : this.y + eps >= p.y;
  }

  le(p) {
    return Math.abs(this.x - p.x) > eps ? this.x <= p.x + eps : this.y <= p.y + eps;
  }

  dot(p) {
    return this.x * p.x + this.y * p.y;
  }

  cross(p) {
    return this.x * p.y - this.y * p.x;
  }

  dist(p = new Point(0, 0)) {
    return Math.hypot(this.x - p.x, this.y - p.y);
  }

  angle() {
    return Math.atan2(this.y, this.x);
  }

  unit() {
    return this.scalaDiv(this.dist());
  }

  rotate(a, p = new Point(0, 0)) {
    const dx = this.x - p.x;
    const dy = this.y - p.y;
    const sin = Math.sin(a);
    const cos = Math.cos(a);
    return new Point(dx * cos - dy * sin + p.x, dx * sin + dy * cos + p.y);
  }

  reflect(line) {
    const u = line.p2.sub(line.p1).unit();
    return line.p1.add(u.scalaMul(this.sub(line.p1).dot(u))).scalaMul(2).sub(this);
  }

  toString() {
    return `(${this.x}, ${this.y})`;
  }
}

class Segment {
  constructor(p1 = new Point(0, 0), p2 = new Point(0, 0)) {
    if (p1.eq(p2)) throw TypeError();
    this.p1 = p1;
    this.p2 = p2;
  }

  len() {
    return this.p1.dist(this.p2);
  }

  dist(p = new Point(0, 0)) {
    const d = this.len() ** 2;
    const t = Math.min(d, Math.max(0, p.sub(this.p1).dot(this.p2.sub(this.p1))));
    return p.sub(this.p1).scalaMul(d).dist(this.p2.sub(this.p1).scalaMul(t)) / d;
  }

  on(p) {
    return this.p1.sub(p).cross(this.p2.sub(p)) === 0 && this.p1.sub(p).dot(this.p2.sub(p)) <= 0;
  }

  intersection(s) {
    const oa = s.p2.sub(s.p1).cross(this.p1.sub(s.p1));
    const ob = s.p2.sub(s.p1).cross(this.p2.sub(s.p1));
    const oc = this.p2.sub(this.p1).cross(s.p1.sub(this.p2));
    const od = this.p2.sub(this.p1).cross(s.p2.sub(this.p1));
    if (oa * ob < 0 && oc * od < 0) return this.p1.scalaMul(ob).sub(this.p2.scalaMul(oa)).scalaDiv(ob - oa);
    const result = [];
    if (s.on(this.p1)) result.push(this.p1);
    if (s.on(this.p2)) result.push(this.p2);
    if (this.on(s.p1)) result.push(s.p1);
    if (this.on(s.p2)) result.push(s.p2);
    if (result.length === 0) return null;
    if (result.length === 1 || result[0].eq(result[1])) return result[0];
    return new Segment(result[0], result[1]);
  }

  move(p) {
    return new Segment(this.p1.add(p), this.p2.add(p));
  }

  rotate(a, p = new Point(0, 0)) {
    return new Segment(this.p1.rotate(a, p), this.p2.rotate(a, p));
  }

  reflect(line) {
    return new Segment(this.p1.reflect(line), this.p2.reflect(line));
  }

  toString() {
    return `${this.p1}~${this.p2}`;
  }
}

class PriorityQueue {
  constructor(f) {
    this.heap1 = [];
    this.heap2 = [];
    this.compare = f;
    this.t = null;
  }

  parentI(i) {
    return Math.floor((i - 1) / 2);
  }

  leftChildI(i) {
    return i * 2 + 1;
  }

  rightChildI(i) {
    return i * 2 + 2;
  }

  childI(i) {
    const leftChildI = this.leftChildI(i);
    const rightChildI = this.rightChildI(i);
    if (leftChildI >= this.heap1.length) return -1;
    if (rightChildI >= this.heap1.length || this.compare(this.heap2[leftChildI], this.heap2[rightChildI]) < 0) return leftChildI;
    return rightChildI;
  }

  push(e1, e2) {
    this.heap1.push(e1);
    this.heap2.push(e2);
    this.heapifyUp();
  }

  pop1() {
    if (this.heap1.length === 0) return null;
    if (this.heap1.length === 1) {
      this.t = this.heap2.pop();
      return this.heap1.pop();
    }
    const top1 = this.heap1[0];
    const bottom1 = this.heap1.pop();
    this.heap1[0] = bottom1;
    const top2 = this.heap2[0];
    const bottom2 = this.heap2.pop();
    this.heap2[0] = bottom2;
    this.heapifyDown();
    this.t = top2;
    return top1;
  }

  pop2() {
    return this.t;
  }

  heapifyUp() {
    let i = this.heap1.length - 1;
    let parentI = this.parentI(i);
    while (parentI !== -1 && this.compare(this.heap2[i], this.heap2[parentI]) < 0) {
      const t1 = this.heap1[parentI];
      this.heap1[parentI] = this.heap1[i];
      this.heap1[i] = t1;
      const t2 = this.heap2[parentI];
      this.heap2[parentI] = this.heap2[i];
      this.heap2[i] = t2;
      i = parentI;
      parentI = this.parentI(i);
    }
  }

  heapifyDown() {
    let i = 0;
    let childI = this.childI(i);
    while (childI !== -1 && this.compare(this.heap2[i], this.heap2[childI]) > 0) {
      const t1 = this.heap1[childI];
      this.heap1[childI] = this.heap1[i];
      this.heap1[i] = t1;
      const t2 = this.heap2[childI];
      this.heap2[childI] = this.heap2[i];
      this.heap2[i] = t2;
      i = childI;
      childI = this.childI(i);
    }
  }

  isEmpty() {
    return this.heap1.length === 0;
  }

  size() {
    return this.heap1.length;
  }

  peek() {
    return this.heap1.length === 0 ? null : [this.heap1[0], this.heap2[0]];
  }
}

function dijkstra(N, E1, E2, sources) {
  const dist = Array(N).fill(Infinity);
  const heap = new PriorityQueue((a, b) => a - b);
  for (const r of sources) {
    dist[r] = 0;
    heap.push(r, 0);
  }

  while (!heap.isEmpty()) {
    const u = heap.pop1();
    const wu = heap.pop2();
    if (wu > dist[u]) continue;
    for (let i = 0; i < 4; i++) {
      const v = E1[i][u];
      const wv = E2[i][u];
      if (isNaN(wv)) continue;
      if (wu + wv >= dist[v]) continue;
      dist[v] = wu + wv;
      heap.push(v, wu + wv);
    }
  }

  return dist;
}

const T = +io.input();

for (let z = 0; z < T; z++) {
  const N = +io.input();
  const R = +io.input();

  const points = new Float64Array(1000000);
  let t = 0;
  const lines = new Uint32Array(N * 4);
  const sources = [];
  const A = Array.from({ length: N }, () => []);

  for (let i = 0; i < N; i++) {
    const sx = +io.input();
    const sy = +io.input();
    const ex = +io.input();
    const ey = +io.input();
    lines[i * 4] = sx;
    lines[i * 4 + 1] = sy;
    lines[i * 4 + 2] = ex;
    lines[i * 4 + 3] = ey;

    const m = +io.input();
    for (let j = 0; j < m; j++) {
      const c = +io.input();
      sources.push(t);
      A[i].push(t);
      points[t * 2] = sx * (1 - c) + ex * c;
      points[t * 2 + 1] = sy * (1 - c) + ey * c;
      t++;
    }
  }

  for (let i = 0; i < N; i++) for (let j = 0; j < i; j++) {
    const sx1 = lines[i * 4];
    const sy1 = lines[i * 4 + 1];
    const ex1 = lines[i * 4 + 2];
    const ey1 = lines[i * 4 + 3];
    const line1 = new Segment(new Point(sx1, sy1), new Point(ex1, ey1));
    const sx2 = lines[j * 4];
    const sy2 = lines[j * 4 + 1];
    const ex2 = lines[j * 4 + 2];
    const ey2 = lines[j * 4 + 3];
    const line2 = new Segment(new Point(sx2, sy2), new Point(ex2, ey2));
    const p = line1.intersection(line2);
    if (p === null) continue;
    A[i].push(t);
    A[j].push(t);
    points[t * 2] = p.x
    points[t * 2 + 1] = p.y;
    t++;
  }

  const Q = +io.input();

  const B = Array.from({ length: Q }, () => []);
  const C = new Uint32Array(Q * 2);
  for (let i = 0; i < Q; i++) {
    const x = +io.input();
    const y = +io.input();
    if (sources.some(e => (x - points[e * 2]) ** 2 + (y - points[e * 2 + 1]) ** 2 <= R ** 2 + eps)) {
      C[i * 2] = 4294967295;
      continue;
    }
    const p = new Point(x, y);
    C[i * 2] = x;
    C[i * 2 + 1] = y;
    if (R === 0) {
      for (let j = 0; j < N; j++) {
        const sx = lines[j * 4];
        const sy = lines[j * 4 + 1];
        const ex = lines[j * 4 + 2];
        const ey = lines[j * 4 + 3];
        const line = new Segment(new Point(sx, sy), new Point(ex, ey));
        if (!line.on(p)) continue;
        A[j].push(t);
        B[i].push(t);
        points[t * 2] = p.x
        points[t * 2 + 1] = p.y;
        t++;
      }
    } else {
      for (let j = 0; j < N; j++) {
        const sx = lines[j * 4];
        const sy = lines[j * 4 + 1];
        const ex = lines[j * 4 + 2];
        const ey = lines[j * 4 + 3];
        const line = new Segment(new Point(sx, sy), new Point(ex, ey));
        const x1 = sx - x;
        const x2 = ex - x;
        const y1 = sy - y;
        const y2 = ey - y;
        const dx = x2 - x1;
        const dy = y2 - y1;
        const drSq = dx ** 2 + dy ** 2;
        const D = x1 * y2 - x2 * y1;
        const dis = R ** 2 * drSq - D ** 2;
        if (dis < 0) continue;
        const s1 = new Point((D * dy + Math.sign(dy + eps) * dx * dis ** 0.5) / drSq + x, (-D * dx + Math.abs(dy) * dis ** 0.5) / drSq + y);
        if ((line.p1.ge(s1) && line.p2.le(s1)) || (line.p2.ge(s1) && line.p1.le(s1))) {
          A[j].push(t);
          B[i].push(t);
          points[t * 2] = s1.x
          points[t * 2 + 1] = s1.y;
          t++;
        }
        const s2 = new Point((D * dy - Math.sign(dy + eps) * dx * dis ** 0.5) / drSq + x, (-D * dx - Math.abs(dy) * dis ** 0.5) / drSq + y);
        if ((line.p1.ge(s2) && line.p2.le(s2)) || (line.p2.ge(s2) && line.p1.le(s2))) {
          A[j].push(t);
          B[i].push(t);
          points[t * 2] = s2.x
          points[t * 2 + 1] = s2.y;
          t++;
        }
      }
    }
  }

  if (t * 2 > 1000000) assert();

  const E1 = Array.from({ length: 4 }, () => new Uint32Array(t).fill(0));
  const E2 = Array.from({ length: 4 }, () => new Float64Array(t).fill(NaN));
  for (let i = 0; i < N; i++) {
    A[i].sort((a, b) => new Point(points[a * 2], points[a * 2 + 1]).ge(new Point(points[b * 2], points[b * 2 + 1])) ? 1 : -1);
    for (let j = 0; j < A[i].length - 1; j++) {
      const u = A[i][j];
      const v = A[i][j + 1];
      const c = Math.hypot(points[u * 2] - points[v * 2], points[u * 2 + 1] - points[v * 2 + 1]);
      const j1 = E2.findIndex(e => isNaN(e[u]));
      E1[j1][u] = v;
      E2[j1][u] = c;
      const j2 = E2.findIndex(e => isNaN(e[v]));
      E1[j2][v] = u;
      E2[j2][v] = c;
    }
  }

  const dijk = dijkstra(t, E1, E2, sources);

  for (let i = 0; i < Q; i++) {
    if (C[i * 2] === 4294967295) {
      io.print(0);
    } else {
      let min = Infinity;
      for (const u of B[i]) min = Math.min(min, dijk[u]);
      io.print(min === Infinity ? -1 : min);
    }
  }
}