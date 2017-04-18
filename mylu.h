#ifndef __MYLU_H
#define __MYLU_H

#include<vector>

#include"base.h"

#define ROW 0
#define COL 1

#define UNVISITED 0
#define VISITED 1

void FormatMatrixIntoCSR(Matrix2D<double> matA,int &n,int &nnz,int **ap,int **ai,double **ax);

void CSR2CSC(int n,int nnz,int *ap,int *ai,double *ax);

void CSC2CSR(int n,int nnz,int *ap,int *ai,double *ax);

void MatrixPermute(int n,int nnz,int *ap,int *ai,double *ax,int *perm);

int DFS(int n,int col,int *luap,int *ludiag,std::vector<int> luai,char *line);

int SymbolPrediction(int n,int k,int *cap,int *cai,int *luap,int *ludiag,std::vector<int> luai,char *nonzero);

int SymbolFactor(int n,int *ap,int *ai,int &lunnz,int *luap,int *ludiag,int **luai);

int PreAnalysis(int n,int nnz,int *ap,int *ai,double *ax,int &lunnz,int *luap,int *ludiag,int **luai,\
		int *rowPerm,int *rowPermInv,int *colPerm,int *colPermInv);

void GPLUFactorize(int n,int *ap,int *ai,double *ax,int *luap,int *ludiag,int *luai,double **luax);

void Lusolve(int n,int *luap,int *ludiag,int *luai,double *luax,double *b,int *rowPerm,int *rowPermInv,int *colPerm,int *colPermInv);

#endif
