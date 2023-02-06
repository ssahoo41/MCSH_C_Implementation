#include <stdio.h>
#include <stdlib.h>

int main() {
    // Write C code here
    int **bin_i = (int*) malloc(sizeof(int*)*37);
    for (int i = 0; i < 37; i ++){
        bin_i[i] = (int*) malloc(sizeof(int)*4);
    }
    printf("Value of bin_i[0][0] %d",bin_i[0][0]);
    for (int i = 0; i < 37; i ++){
        free(bin_i[i]);
    }
    free(bin_i);
    return 0;
}