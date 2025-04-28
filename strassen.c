#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void add_matrix(double* C, double* A, double* B, int n) {
    for (int i = 0; i < n * n; i++) {
        C[i] = A[i] + B[i];
    }
}

void sub_matrix(double* C, double* A, double* B, int n) {
    for (int i = 0; i < n * n; i++) {
        C[i] = A[i] - B[i];
    }
}

void strassen(double* C, double* A, double* B, int n) {
    int new_n = n / 2;
    int size = new_n * new_n;

    double *A00 = A;
    double *A01 = A + new_n;
    double *A10 = A + n*new_n;
    double *A11 = A + n*new_n + new_n;
    double *B00 = B;
    double *B01 = B + new_n;
    double *B10 = B + n*new_n;
    double *B11 = B + n*new_n + new_n;
    double *C00 = C;
    double *C01 = C + new_n;
    double *C10 = C + n*new_n;
    double *C11 = C + n*new_n + new_n;

    double *M0 = (double*)malloc(sizeof(double)*size);
    double *M1 = (double*)malloc(sizeof(double)*size);
    double *M2 = (double*)malloc(sizeof(double)*size);
    double *M3 = (double*)malloc(sizeof(double)*size);
    double *M4 = (double*)malloc(sizeof(double)*size);
    double *M5 = (double*)malloc(sizeof(double)*size);
    double *M6 = (double*)malloc(sizeof(double)*size);
    double *T1 = (double*)malloc(sizeof(double)*size);
    double *T2 = (double*)malloc(sizeof(double)*size);

    // M0 = (A00 + A11)*(B00 + B11)
    add_matrix(T1, A00, A11, new_n);
    add_matrix(T2, B00, B11, new_n);
    strassen(M0, T1, T2, new_n);

    // M1 = (A10 + A11)*B00
    add_matrix(T1, A10, A11, new_n);
    strassen(M1, T1, B00, new_n);

    // M2 = A00*(B01 - B11)
    sub_matrix(T1, B01, B11, new_n);
    strassen(M2, A00, T1, new_n);

    // M3 = A11*(B10 - B00)
    sub_matrix(T1, B10, B00, new_n);
    strassen(M3, A11, T1, new_n);

    // M4 = (A00 + A01)*B11
    add_matrix(T1, A00, A01, new_n);
    strassen(M4, T1, B11, new_n);

    // M5 = (A10 - A00)*(B00 + B01)
    sub_matrix(T1, A10, A00, new_n);
    add_matrix(T2, B00, B01, new_n);
    strassen(M5, T1, T2, new_n);

    // M6 = (A01 - A11)*(B10 + B11)
    sub_matrix(T1, A01, A11, new_n);
    add_matrix(T2, B10, B11, new_n);
    strassen(M6, T1, T2, new_n);

    
    for (int i=0; i < new_n; i++) {
        for (int j=0; j <new_n; j++) {
            C00[i*n+j]= M0[i*new_n+j]+M3[i*new_n+j]-M4[i*new_n+j] + M6[i*new_n+j];
            C01[i*n+j]= M2[i*new_n+j] +M4[i*new_n+j];
            C10[i*n+j] = M1[i*new_n+j]+ M3[i*new_n+j];
            C11[i*n+j]= M0[i*new_n+j] - M1[i*new_n+j] + M2[i*new_n+j] + M5[i*new_n+j];
        }
    }

    free(M0); 
    free(M1); 
    free(M2); 
    free(M3);
    free(M4); 
    free(M5); 
    free(M6);
    free(T1); 
    free(T2);
}
