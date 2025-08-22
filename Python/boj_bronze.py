Y, M, D = map(int, input().split('-'))
N = int(input())
D += N
if D > 30:
    M += 1
    D -= 30
if M > 12:
    Y += 1
    M -= 12
print(f"{Y}-{M:02}-{D:02}")
