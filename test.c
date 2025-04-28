#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

double get_time_in_seconds() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}


void strassen(double* C, double* A, double* B, int n);

void normal_gemm(double* C, double* A, double* B, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i*n + k] * B[k*n + j];
            }
            C[i*n + j] = sum;
        }
    }
}


void random_matrix(double* M, int n) {
    for (int i = 0; i < n * n; i++) {
        M[i] = (double)(rand() % 10);
    }
}

void print_matrix(double* M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%5.1f ", M[i*n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main() {
    

    FILE *fp = fopen("output.txt", "w");
    fprintf(fp, "n, normal_time_sec, strassen_time_sec, max_difference\n");

    for (int dim = 64; dim < 2048; dim+=64) {
        int n = dim;

        double* A = (double*)malloc(sizeof(double) * n * n);
        double* B = (double*)malloc(sizeof(double) * n * n);
        double* C1 = (double*)malloc(sizeof(double) * n * n);
        double* C2 = (double*)malloc(sizeof(double) * n * n);

        
        random_matrix(A, n);
        random_matrix(B, n);

        double start, end, normal_time, strassen_time;

        //gemm
        start = get_time_in_seconds();
        normal_gemm(C1, A, B, n);
        end = get_time_in_seconds();
        normal_time = end - start;

        //strassen
        start = get_time_in_seconds();
        strassen(C2, A, B, n);
        end = get_time_in_seconds();
        strassen_time = end - start;

    
        double max_diff = 0.0;
        for (int i = 0; i < n * n; i++) {
            double diff = fabs(C1[i] - C2[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        fprintf(fp, "%d, %.6f, %.6f, %e\n", n, normal_time, strassen_time, max_diff);

        free(A);
        free(B);
        free(C1);
        free(C2);
    }

    fclose(fp);

    return 0;
}
