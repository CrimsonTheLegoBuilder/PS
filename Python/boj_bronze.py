S = input()
A, B = 0, 0
f = 1
for i in range(0, len(S), 2):
    if S[i] == 'A':
        A += int(S[i + 1])
    if S[i] == 'B':
        B += int(S[i + 1])
    if A == B == 10:
        f = 0
    if f:
        if A == 11:
            print("A")
            break
        if B == 11:
            print("B")
            break
    if abs(A - B) > 1:
        print("A" if A > B else "B")
        break
