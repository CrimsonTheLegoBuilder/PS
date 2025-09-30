//Struct Pos { int x, y, i }
//Struct Pdd - 보로노이 다이어그램 제작에만 쓰임
//Class or Struct voronoi diagram - Seed, Cell :: 재홍이가 접근할 일 없도록 구성
//Class or Struct KDtree - Tree, Pos :: 재홍이가 접근할 일 없도록 구성
//입력 및 전처리는 전부 기하모듈에서 일어나야하는데 후에 하는 접근은 전부 다른 자료구조들에서 행해짐
//전처리를 전부 끝내놓고 나면 기하모듈에 접근할 일 없도록 해야함
//Class or Struct Query - int t, Pos s, Pos e, int u, int p :: 타입, 점 2개, 부모별 번호. 순례끝점 위치는 KDtree 조회 후 i 에 번호 기록