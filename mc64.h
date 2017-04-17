#ifndef __MC64_H
#define __MC64_H

#define INF 1.0E50

#define ABS(x) ((x<0)?-(x):(x))

class MinHeap{
private:
	void Reheap(double *d);
public:
	MinHeap(long size=0);
	~MinHeap();
	void Add(long h, double *d);
	long PopRoot(double *d);
	void Clear();
	void Resize(long size);
public:
	long max;
	long n;	//number of elements in heap
	long *data;
};

long MaxMatch(long n, long *ap, long *ai, long *perm, long *permInv);

long MC64(long n, long nnz, long *ap, long *ai, double *ax, long *perm, long *permInv);

#endif
