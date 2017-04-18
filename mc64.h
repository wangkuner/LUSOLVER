#ifndef __MC64_H
#define __MC64_H

#define INF 1.0E50

#define ABS(x) ((x<0)?-(x):(x))

class MinHeap{
private:
	void Reheap(double *d);
public:
	MinHeap(int size=0);
	~MinHeap();
	void Add(int h, double *d);
	int PopRoot(double *d);
	void Clear();
	void Resize(int size);
public:
	int max;
	int n;	//number of elements in heap
	int *data;
};

int MaxMatch(int n, int *ap, int *ai, int *perm, int *permInv);

int MC64(int n, int nnz, int *ap, int *ai, double *ax, int *perm, int *permInv);

void MC21(int n, int *ap, int *ai, int *perm, int *permInv);

#endif
