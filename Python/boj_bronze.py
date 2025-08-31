for i in range(2992, 10000):
    a = i
    b = int(str(i), 12)
    c = int(str(i), 16)
    # print(a, b, c)
    A = 0
    while a:
        # print(a)
        A += a % 10
        a = a // 10
    B = 0
    while b:
        # print(b)
        B += b % 12
        b = b // 12
    C = 0
    while c:
        # print(c)
        C += c % 16
        c = c // 16
    if A == B == C:
        print(i)