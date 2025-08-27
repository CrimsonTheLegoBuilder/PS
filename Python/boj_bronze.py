N, K = map(int, input().split())
L = []
for i in range(1, K + 1):
    n = int(str(N * i)[::-1])
    L.append(n)
L.sort()
print(L[-1])

