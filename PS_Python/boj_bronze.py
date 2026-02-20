import sys

N = int(sys.stdin.readline())
MOD = 1000000

def mul(a, b):
    c = [[0, 0], [0, 0]]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                c[i][j] += a[i][k] * b[k][j]
            c[i][j] %= MOD
    return c

def power(a, n):
    res = [[1, 0], [0, 1]]
    while n > 0:
        if n % 2 == 1:
            res = mul(res, a)
        a = mul(a, a)
        n //= 2
    return res

if N == 0:
    print(0)
else:
    M = [[1, 1], [1, 0]]
    ret = power(M, N)

    print(ret[0][1])