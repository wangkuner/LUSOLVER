#ifndef __MYLU_H
#define __MYLU_H

#include<vector>

#include"base.h"

#define ROW 0
#define COL 1

#define UNVISITED 0
#define VISITED 1

void FormatMatrixIntoCSR(Matrix2D<double> matA,long &n,long &nnz,long **ap,long **ai,double **ax);

void CSR2CSC(long n,long nnz,long *ap,long *ai,double *ax);

void CSC2CSR(long n,long nnz,long *ap,long *ai,double *ax);

void MatrixPermute(long n,long nnz,long *ap,long *ai,double *ax,long *perm);

int DFS(long n,long col,long *lucap,long *ludiag,std::vector<long> lucai,char *line);

int SymbolPrediction(long n,long k,long *cap,long *cai,long *lucap,long *ludiag,std::vector<long> lucai,char *nonzero);

int SymbolFactor(long n,long nnz,long *cap,long *cai,double *cax,long &lunnz,long *luap,long *ludiag,long **luai,long *rowPerm,long *rowPermInv,long *colPerm,long *colPermInv);

void GPLUFactorize(long n,long &nnz,long *cap,long *cdiag,long **cai,double **cax);

void MyNicsluSolve(long n,long *cap,long *cdiag,long *cai,double *cax,double *b,long *rowPerm,long *rowPermInv,long *colPerm,long *colPermInv);

#endif
