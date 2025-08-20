A, B, C = map(int, input().split())
assert A <= B <= C
D = (A + B + C) // 3
assert D * 3 == A + B + C
T = 0
if D == B:
    T = (C - D) * 2
elif D < B:
    T = B - D
    T += (C - D) * 2
elif D > B:
    T = D - B
    C -= T
    T += (C - D) * 2
print(T)