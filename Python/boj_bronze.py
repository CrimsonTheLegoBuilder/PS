N = int(input())
for i in range(N):
    print(f"Case {i + 1}: ", end='')
    S = input()
    parts = S.split()
    a = int(parts[0])
    b = int(parts[2])
    c = int(parts[4])
    if '+' in S:
        d = a + b
        if d == c:
            print("YES")
        else:
            print("NO")
    else:
        d = a - b
        if d == c:
            print("YES")
        else:
            print("NO")