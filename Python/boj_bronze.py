while 1:
    try:
        N = list(map(int, input().split()))
        for i in range(len(N)):
            k = N[i]
            if i > 0:
                k *= N[i - 1]
            if i < len(N) - 1:
                k *= N[i + 1]
            print(k, end=' ')
        # print()
    except EOFError:
        break
