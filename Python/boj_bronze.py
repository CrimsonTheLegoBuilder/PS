for _ in range(int(input())):
    p = input()
    q = input()
    t = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            t += 1
    print("Hamming distance is", t, end='')
    print(".")