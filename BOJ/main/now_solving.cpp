#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>

int main() {
    int stack[100000] = { 0, };
    int array[100000] = { 0, };
    char act[100000];
    int N, top = -1, element = 1;
    scanf("%d", &N);
    for (int i = 0; i < N; i++) {
        scanf("%d", &array[i]);
    }
    int i = 0;
    int tmp = 0;
    while (1) {
        if (i < N) {
            if (array[i] < element) {
                if (stack[top] != array[i]) {
                    printf("NO\n");
                    break;
                }
            }
            while (array[i] >= element) {
                stack[++top] = element++;
                act[tmp++] = '+';

            }
            if (stack[top] == array[i]) {
                i++;
                act[tmp++] = '-';

                stack[top--] = 0;
            }
            else if (stack[top] != array[i]) {
                if (top == -1) break;
                if (i >= N) {
                    printf("NO\n");
                    break;
                }
            }
        }
        else {
            if (top == -1) {
                for (int n = 0; n < tmp; n++) {
                    printf("%c\n", act[n]);
                }
                break;
            }
            else {
                printf("NO\n");
                break;
            }
        }
    }
    return 0;
}