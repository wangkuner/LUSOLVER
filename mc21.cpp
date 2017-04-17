#include<iostream>
#include<string>

using namespace std;

void MC21(int n, int *ap, int *ai, int *perm, int *permInv);

void MC21(int n, int *ap, int *ai, int *perm,int *permInv)
{
	int *p = new int[2 * n]{ -1 };
	int *amask = new int[n] {0};
	size_t size{ n*sizeof(int) };

	memset(perm, -1, size);
	memset(permInv, -1, size);

	int col{ 0 };
	int row{ 0 };
	int top{ 0 };
	int iap{ -1 };
	int flag{ 0 };

	for (int i = 0; i < n; i++)
	{
		size = n*sizeof(int);
		memset(amask, 0, size);
		size = 2 * n*sizeof(int);
		memset(p, -1, size);

		col = i;
		p[col] = -1;
		row = ap[col];
		iap = -1;
		flag = 0;
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
					flag = 1;
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
		}
		else
		{
			cout << "no permutation matrix" << endl;
			break;
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
}