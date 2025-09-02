for _ in range(int(input())):
    S = list(map(str, input().split()))
    print("god", end='')
    for x in S[1:]:
        print(x, end='')
    print()