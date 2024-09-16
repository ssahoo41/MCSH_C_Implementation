#include <stdlib.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <math.h>

int main(){

    double *params_di = (double*) malloc(sizeof(double)*4);
    FILE *in_file = fopen("GMP.txt", "r");

    if (NULL == in_file) printf("File cannont be opened\n");

    for (int i =0; i < 4; i++){
        fscanf(in_file,"%lf",&params_di[i]);
        printf("%f,",params_di[i]);
    }
    printf("\n");
    double sigmas[3] = {params_di[0],params_di[1],params_di[2]};
    double cutoff = params_di[3];

    double **params_d = (double**) malloc(sizeof(double*)*15);
    for (int i = 0; i < 15; i++) params_d[i] = (double*) malloc(sizeof(double)*6);

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 3; j++){
            params_d[i*j+j][0] = sigmas[j];
            params_d[i*j+j][1] = 1.0;
            params_d[i*j+j][2] = pow(1/(sigmas[j]*sqrt(2.0*M_PI)),3);
            params_d[i*j+j][3] = 1/(2*pow(sigmas[j],2));
            params_d[i*j+j][4] = cutoff;
            params_d[i*j+j][5] = 1.0;
            printf("%f,%f,%f,%f,%f,%f\n",params_d[i*j+j][0],params_d[i*j+j][1],params_d[i*j+j][2],params_d[i*j+j][3],params_d[i*j+j][4],params_d[i*j+j][5]);

        }
    }

    return 0;
}

/*

int main() {
    // Array of pointers
    int *params_ii[3];
    // Single pointer to array with memory allocated for 3 integers
    int *params_ip = (int*) malloc(sizeof(int)*3);

    for (int i =0; i < 3; i++){
        printf("%d,",params_ii[i]);
    }
    printf("\n");
    for (int i =0; i < 3; i++){
        printf("%d,",params_ip[i]);
    }
    return 0;
}
*/