#include<iostream>
#include<iomanip>
#include<cstring>
#include<cmath>
#include"mc64.h"

using namespace std;

MinHeap::MinHeap(int size)
{
	max = size;
	n = 0;

	if ((max > 0))
	{
		data = new int[max] {0};
	}
	else
	{
		data = nullptr;
	}
}

MinHeap::~MinHeap()
{
	max = 0;
	n = 0;

	if (data != nullptr)
	{
		delete[] data;
		data = nullptr;
	}
}

void MinHeap::Reheap(double *d)
{
	int rowRoot = data[0];

	int rootIndex{ 0 };
	int smallChildIndex{ 0 };
	int leftChildIndex{ 0 };
	int rightChildIndex{ 0 };

	bool done{ false };

	leftChildIndex = 2 * (rootIndex + 1) - 1;
	rightChildIndex = leftChildIndex + 1;
	smallChildIndex = leftChildIndex;
	while (!done&&smallChildIndex < n)
	{
		if (rightChildIndex < n&&d[data[rightChildIndex]] < d[data[smallChildIndex]])
		{
			smallChildIndex = rightChildIndex;
		}

		if (d[rowRoot] > d[data[smallChildIndex]])
		{
			data[rootIndex] = data[smallChildIndex];
			rootIndex = smallChildIndex;
			leftChildIndex = 2 * (rootIndex + 1) - 1;
			rightChildIndex = leftChildIndex + 1;
			smallChildIndex = leftChildIndex;
		}
		else
		{
			done = true;
		}
	}

	data[rootIndex] = rowRoot;
}

void MinHeap::Add(int h, double *d)
{
	int index{ n };
	int parent{ (index + 1) / 2 - 1 };
	
	n++;
	while (parent >= 0 && (d[h] < d[data[parent]]))
	{
		data[index] = data[parent];
		index = parent;
		parent = (index + 1) / 2 - 1;
	}

	data[index] = h;
}

int MinHeap::PopRoot(double *d)
{
	int row{ data[0] };

	if (n > 0)
	{
		data[0] = data[n - 1];
		n--;
		Reheap(d);
		return row;
	}

	return -1;
}

void MinHeap::Clear()
{
	n = 0;
}

void MinHeap::Resize(int size)
{
	if (size != max)
	{
		max = size;
		n = 0;
		if (data != nullptr)
		{
			delete[] data;
			data = nullptr;
		}
		if (max > 0)
		{
			data = new int[max] {0};
		}
	}
}

int MaxMatch(int n, int *ap, int *ai, int *perm, int *permInv)
{
	int *p = new int[2 * n]{ -1 };
	char *amask = new char[n] {0};
	size_t size{ n*sizeof(int) };

	for (int i = 0; i < n; i++)
	{
		*(perm + i) = -1;
		*(permInv + i) = -1;
	}

	int col{ 0 };
	int row{ 0 };
	int top{ 0 };
	int iap{ -1 };
	bool flag{ false };
	int finish{ 0 };

	for (int i = 0; i < n; i++)
	{
		size = n*sizeof(char);
		memset(amask, 0, size);
		//size = 2 * n*sizeof(int);
		//memset(p, 0, size);

		col = i;
		p[col] = -1;
		row = ap[col];
		iap = -1;
		flag = false;
		top = 0;

		while (row < ap[col + 1] || top>0)	//find an augmenting path form col i
		{
			if (row < ap[col + 1])
			{
				if (permInv[ai[row]] < 0)	//row j unmatched
				{
					iap = ai[row];
					break;
				}
				else if (!amask[ai[row]])	//row j unvisited
				{
					amask[ai[row]] = 1;
					top++;
					p[2 * top] = col;
					p[2 * top + 1] = row;
					col = permInv[ai[row]];
					row = ap[col];
					flag = true;
				}
				else
				{
					row++;
				}
			}
			else if (top > 0)
			{
				col = p[2 * top];
				row = p[2 * top + 1] + 1;
				top--;
			}
		}

		if (iap >= 0)
		{
			if (flag)
			{
				while (top > 0)
				{
					perm[p[2 * top]] = ai[p[2 * top + 1]];
					permInv[ai[p[2 * top + 1]]] = p[2 * top];
					top--;
				}
				perm[col] = iap;
				permInv[iap] = col;
			}
			else
			{
				perm[i] = iap;
				permInv[iap] = i;
			}
			finish++;
		}
		else
		{
			continue;
		}
	}

	if (p != nullptr)
	{
		delete[] p;
		p = nullptr;
	}
	if (amask != nullptr)
	{
		delete[] amask;
		amask = nullptr;
	}
	return finish;
}

int MC64(int n, int nnz, int *ap, int *ai, double *ax, int *perm, int *permInv)
{
	double *cx{ nullptr };		//nnz
	double *cxp{ nullptr };		//nnz,cx prime for update
	double *amax{ nullptr };	//n,maxminium element of column
	double *u{ nullptr };		//n
	double *v{ nullptr };		//n
	int *bset{ nullptr };		//n
	int *prevCol{ nullptr };	//n
	MinHeap qset(n);
	int *qmask{ nullptr };		//n

	double *d{ nullptr };		//n,length of row i to column jo

	int *apLocal{ nullptr };	//n+1
	int *aiLocal{ nullptr };	//not bigger than nnz

	int maxMatch{ 0 };
	int col{ 0 };
	double lsp{ 0.0 };			//length of shortest path 
	double lsap{ 0.0 };			//length of shortest augmenting path
	int rsap{ -1 };				//exist augmenting path from row isap to root column
	int csap{ -1 };
	bool flag{ true };
	double dnew{ 0.0 };
	int row{ 0 };
	int tmpRow{ 0 };
	int tmpCol{ 0 };

	cx = new double[2 * nnz + 4 * n]{ 0.0 };
	if (cx == nullptr)
	{
		cout << "in file \"" << __FILE__ << "\" line " << __LINE__ << ":" << "allocating memory failed!" << endl;
		return -2;
	}
	cxp = cx + nnz;
	amax = cxp + nnz;
	u = amax + n;
	v = u + n;
	d = v + n;

	bset = new int[3 * n + nnz + n + 1]{ 0 };
	if (bset == nullptr)
	{
		cout <<"in file \""<<__FILE__<<"\" line "<<__LINE__<<":"<< "allocating memory failed!" << endl;
		return -2;
	}
	prevCol = bset + n;
	qmask = prevCol + n;
	apLocal = qmask + n;
	aiLocal = apLocal + n + 1;
	//reset the match
	for (int i = 0; i < n; i++)
	{
		*(perm + i) = -1;
		*(permInv + i) = -1;
	}

	//find maxminium element of every column
	for (int i = 0; i < n; i++)
	{
		amax[i] = ABS(ax[ap[i]]);
		for (int j = ap[i]; j < ap[i + 1]; j++)
		{
			if (amax[i] < ABS(ax[j]))
			{
				amax[i] = ABS(ax[j]);
			}
		}
	}
	//reduced matrix
	apLocal[0] = 0;
	for (int i = 0; i < n; i++)
	{
		apLocal[i + 1] = apLocal[i];
		for (int j = ap[i]; j < ap[i + 1]; j++)
		{
			cx[j] = log(amax[i]) - log(ABS(ax[j]));
			cxp[j] = cx[j];
			if (cx[j] < 1.0e-32)
			{
				aiLocal[apLocal[i + 1]] = ai[j];
				apLocal[i + 1]++;
			}
		}
	}
	//fine initial max match
	maxMatch = MaxMatch(n, apLocal, aiLocal, perm, permInv);

	cout << setw(10) << "perm:";
	for (int i = 0; i < n; i++)
	{
		cout << setw(10) << perm[i];
	}
	cout << endl;
	cout << setw(10) << "permInv:";
	for (int i = 0; i < n; i++)
	{
		cout << setw(10) << permInv[i];
	}
	cout << endl;

	if (maxMatch < n)
	{
		for (int jo = 0; jo < n; jo++)
		{
			if (perm[jo] < 0)	//column j0 unmatched
			{
				memset(bset, 0, n*sizeof(int));
				memset(qmask, 0, n*sizeof(int));
				qset.Clear();
				for (int i = 0; i < n; i++)
				{
					d[i] = INF;
				}
				lsp = 0.0;
				lsap = INF;
				rsap = -1;
				col = jo;
				csap = col;
				prevCol[col] = -1;

				flag = true;
				while (flag)
				{
					for (int i = ap[col]; i < ap[col + 1]; i++)
					{
						row = ai[i];
						if (!bset[row])
						{
							dnew = lsp + cxp[i];
							if (dnew < lsap)
							{
								if (permInv[row] < 0)
								{
									lsap = dnew;
									rsap = row;
									csap = col;
								}
								else
								{
									if (dnew < d[row])
									{
										d[row] = dnew;
										prevCol[permInv[row]] = col;
										if (!qmask[row])
										{
											qset.Add(row, d);
											qmask[row] = 1;
										}
									}
								}
							}
						}
					}
					if (qset.n < 1)
					{
						break;
					}
					row = qset.PopRoot(d);
					lsp = d[row];
					if (lsap < lsp)
					{
						break;
					}
					bset[row] = 1;
					col = permInv[row];
				}
				//next augmenting path
				if (rsap >= 0)
				{
					//add augmenting path into match M
					tmpRow = perm[csap];
					//tmpCol = col;
					permInv[rsap] = csap;
					perm[csap] = rsap;
					col = csap; 
					col = prevCol[col];
					while (col >= 0)
					{
						tmpCol = perm[col];
						perm[col] = tmpRow;
						permInv[tmpRow] = col;
						col = prevCol[col];
						tmpRow = tmpCol;
					}
					maxMatch++;
					if (maxMatch >= n)
					{
						break;
					}
					
					//update u
					for (int i = 0; i < n; i++)
					{
						if (bset[i])
						{
							u[i] += d[i] - lsap;
						}
					}
					//update v
					for (int i = 0; i < n; i++)
					{
						for (int j = ap[i]; j < ap[i + 1]; j++)
						{
							row = ai[j];
							if (row == perm[i])
							{
								v[i] = cx[j] - u[row];
							}
						}
					}
					//update C
					for (int i = 0; i < n; i++)
					{
						for (int j = ap[i]; j < ap[i + 1]; j++)
						{
							row = ai[j];
							cxp[j] = cx[j] - u[row] - v[i];
						}
					}
				}
				else
				{
					break;
				}
			}
			col = jo;
		}
		if (maxMatch >= n)
		{
			col = -1;
		}
	}
	else
	{
		col = -1;
	}

	if (cx != nullptr)
	{
		delete[] cx;
		cx = nullptr;
	}

	if (bset != nullptr)
	{
		delete[] bset;
		bset = nullptr;
	}

	return col;
}
