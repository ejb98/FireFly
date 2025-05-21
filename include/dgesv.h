#ifndef DGESV_H
#define DGESV_H

extern void dgesv_(int* n, int* nrhs, double* a, int* lda,
                   int* ipiv, double* b, int* ldb, int* info);

#endif