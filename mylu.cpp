#include<iostream>
#include<cstring>
#include<iterator>

#include"mylu.h"
#include"mc64.h"
#include"amd.h"

using std::cout;
using std::endl;
using std::vector;



void FormatMatrixIntoCSR(Matrix2D<double> matA,int &n,int &nnz,int **ap,int **ai,double **ax)
{
	unsigned int num{0};

	if(matA.row!=matA.col)
	{
		cout<<"matrix is not square"<<endl;
		exit(0);
	}

	n=matA.row;
	
	if(*ap!=nullptr)
	{
		delete[] *ap;
	}
	*ap=new int[n+1]{0};

	for(unsigned int i=0;i<matA.row;i++)
	{
		for(unsigned int j=0;j<matA.col;j++)
		{
			if(!(*(matA.Mat+i*matA.col+j)<1e-6&&*(matA.Mat+i*matA.col+j)>-1e-6))
			{
				num++;
			}
		}
		*(*ap+i+1)=num;
	}

	nnz=num;
	if(*ax!=nullptr)
	{
		delete[] *ax;
	}
	*ax=new double[nnz]{0};
	if(*ai!=nullptr)
	{
		delete[] *ai;
	}
	*ai=new int[nnz]{0};

	num=0;

	for(unsigned int i=0;i<matA.row;i++)
	{
		for(unsigned int j=0;j<matA.col;j++)
		{
			if(!(*(matA.Mat+i*matA.col+j)<1e-6&&*(matA.Mat+i*matA.col+j)>-1e-6))
			{
				*(*ax+num)=*(matA.Mat+i*matA.col+j);
				*(*ai+num)=j;
				num++;
			}
		}
	}
}

void CSR2CSC(int n,int nnz,int *ap,int *ai,double *ax)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配CSC存储空间"<<endl;
		exit(0);
	}

	size_t size=0;

	int *nonzeropercol=new int[n]{0};

	int *cap=new int[n+1]{0};
	int *cai=new int[nnz]{0};
	double *cax=new double[nnz]{0};

	for(int i=0;i<n;i++)
	{
		for(int j=ap[i];j<ap[i+1];j++)
		{
			cap[ai[j]+1]++;
		}
	}

	for(int i=0;i<n;i++)
	{
		cap[i+1]+=cap[i];
	}

	for(int i=0;i<n;i++)
	{
		for(int j=ap[i];j<ap[i+1];j++)
		{
			cai[cap[ai[j]]+nonzeropercol[ai[j]]]=i;
			cax[cap[ai[j]]+nonzeropercol[ai[j]]]=ax[j];
			nonzeropercol[ai[j]]++;
		}
	}

	size=(n+1)*sizeof(int);
	memcpy(ap,cap,size);
	size=nnz*sizeof(int);
	memcpy(ai,cai,size);
	size=nnz*sizeof(double);
	memcpy(ax,cax,size);

	if(nonzeropercol!=nullptr)
	{
		delete[] nonzeropercol;
		nonzeropercol=nullptr;
	}
	if(cap!=nullptr)
	{
		delete[] cap;
		cap=nullptr;
	}
	if(cai!=nullptr)
	{
		delete[] cai;
		cai=nullptr;
	}
	if(cax!=nullptr)
	{
		delete[] cax;
		cax=nullptr;
	}
}

void CSC2CSR(int n,int nnz,int *ap,int *ai,double *ax)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配CSR存储空间"<<endl;
		exit(0);
	}

	size_t size=0;

	int *nonzeropercol=new int[n]{0};

	int *cap=new int[n+1]{0};
	int *cai=new int[nnz]{0};
	double *cax=new double[nnz]{0};

	size=(n+1)*sizeof(int);
	memcpy(cap,ap,size);
	memset(ap,0,size);
	size=nnz*sizeof(int);
	memcpy(cai,ai,size);
	size=nnz*sizeof(double);
	memcpy(cax,ax,size);

	int *nonzeroperrow=new int[n]{0};

	for(int i=0;i<n;i++)
	{
		for(int j=cap[i];j<cap[i+1];j++)
		{
			ap[cai[j]+1]++;
		}
	}

	for(int i=0;i<n;i++)
	{
		ap[i+1]+=ap[i];
	}

	for(int i=0;i<n;i++)
	{
		for(int j=cap[i];j<cap[i+1];j++)
		{
			ai[ap[cai[j]]+nonzeroperrow[cai[j]]]=i;
			ax[ap[cai[j]]+nonzeroperrow[cai[j]]]=cax[j];
			nonzeroperrow[cai[j]]++;
		}
	}

	if(nonzeroperrow!=nullptr)
	{
		delete[] nonzeroperrow;
		nonzeroperrow=nullptr;
	}
	if(cap!=nullptr)
	{
		delete[] cap;
		cap=nullptr;
	}
	if(cai!=nullptr)
	{
		delete[] cai;
		cai=nullptr;
	}
	if(cax!=nullptr)
	{
		delete[] cax;
		cax=nullptr;
	}
}

void MatrixPermute(int n,int nnz,int *ap,int *ai,double *ax,int *perm)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配存储空间"<<endl;
		exit(0);
	}

	size_t size=0;
	int p=0;

	int *tap=new int[n+1]{0};
	int *tai=new int[nnz]{0};
	double *tax=new double[nnz]{0};

	size=(n+1)*sizeof(int);
	memcpy(tap,ap,size);
	size=(nnz)*sizeof(int);
	memcpy(tai,ai,size);
	size=(nnz)*sizeof(double);
	memcpy(tax,ax,size);

	ap[0]=0;
	for(int i=0;i<n;i++)
	{
		ap[i+1]=ap[i]+tap[perm[i]+1]-tap[perm[i]];
	}

	for(int i=0;i<n;i++)
	{
		p=ap[i];
		for(int j=tap[perm[i]];j<tap[perm[i]+1];j++)
		{
			ai[p]=tai[j];
			ax[p]=tax[j];
			p++;
		}
	}

	if(tap!=nullptr)
	{
		delete[] tap;
		tap=nullptr;
	}
	if(tai!=nullptr)
	{
		delete[] tai;
		tai=nullptr;
	}
	if(tax!=nullptr)
	{
		delete[] tax;
		tax=nullptr;
	}
}

int DFS(int n,int col,int *luap,int *ludiag,std::vector<int> luai,char *line)
{
	/*递归方法速度太慢
	if(line[colki]==VISITED)
	{
		return 0;
	}

	int cstart{0};
	int cend{0};

	line[colki]=VISITED;

	if(luap[colki+1]-ludiag[colki]>1)
	{
		cstart=ludiag[colki]+1;
		cend=luap[colki+1];
		for(int i=cstart;i<cend;i++)
		{
			DFS(luai[i],luap,ludiag,luai,line);
		}
	}
	else
	{
		return 0;
	}
	*/
	int coli{ludiag[col]+1};
	vector<int> stack(2*n);
	int top=-1;

	
	line[col]=VISITED;

	while(coli<luap[col+1] || (!(top==-1)))
	{
		if(coli<luap[col+1])
		{
			while(coli<luap[col+1])
			{
				if(line[luai[coli]]!=VISITED)
				{
					stack[++top]=col;
					stack[++top]=coli;
					col=luai[coli];
					coli=ludiag[col]+1;
					line[col]=VISITED;
					break;
				}
				coli++;
			}
		}
		else
		{
			coli=stack[top--];
			//stack.pop_back();
			col=stack[top--];
			//stack.pop_back();
			while(coli<luap[col+1])
			{
				if(line[luai[coli]]!=VISITED)
				{
					//stack.push_back(col);
					//stack.push_back(coli);
					stack[++top]=col;
					stack[++top]=coli;
					col=luai[coli];
					coli=ludiag[col]+1;
					line[col]=VISITED;
					break;
				}
				coli++;
			}
		}
	}

	return 1;
}

int SymbolPrediction(int n,int k,int *ap,int *ai,int *luap,int *ludiag,std::vector<int> luai,char *nonzero)
{
	char *line=new char[n]{0};

	int cstart=ap[k];
	int cend=ap[k+1];
	int ln=0;

	for(int i=cstart;i<cend;i++)
	{
		DFS(n,ai[i],luap,ludiag,luai,line);
		
		for(ln=0;ln<n;ln++)
		{
			if(line[ln]==VISITED)
			{
				nonzero[ln]=VISITED;
				line[ln]=UNVISITED;
			}
		}
	}

	if(line!=nullptr)
	{
		delete[] line;
		line=nullptr;
	}
	return 1;
}

int SymbolFactor(int n,int *ap,int *ai,int &lunnz,int *luap,int *ludiag,int **luai)
{
	vector<int>tmpluai;

	char *nonzero=new char[n]{0};

	int len{0};
	int diag{0};

	if(!nonzero)
	{
		cout<<"file \""<<__FILE__<<"\" line "<<__LINE__<<":"<<"allocating memory failed"<<endl;
		return 0;
	}

	//symbol prediction, get the parttern of LU
	for(int k=0;k<n;k++)
	{
		SymbolPrediction(n,k,ap,ai,luap,ludiag,tmpluai,nonzero);

		len=0;
		diag=0;
		for(int i=0;i<n;i++)
		{
			if(i<=k)
			{
				if(nonzero[i]!=0)
				{
					len++;
					diag++;
					tmpluai.push_back(i);
					nonzero[i]=0;
				}
			}
			else
			{
				if(nonzero[i]!=0)
				{
					len++;
					tmpluai.push_back(i);
					nonzero[i]=0;
				}
			}
		}
		luap[k+1]=luap[k]+len;
		ludiag[k]=luap[k]+diag-1;
	}

	lunnz=tmpluai.size();

	if(*luai!=nullptr)
	{
		delete[] *luai;
	}
	*luai=new int[lunnz];
	if(!(*luai))
	{
		cout<<"file \""<<__FILE__<<"\" line "<<__LINE__<<":"<<"allocating memory failed"<<endl;
		return 0;
	}

	for(int i=0;i<lunnz;i++)
	{
		(*luai)[i]=tmpluai[i];
	}

	if(nonzero!=nullptr)
	{
		delete[] nonzero;
		nonzero=nullptr;
	}

	return 1;
}

int PreAnalysis(int n,int nnz,int *ap,int *ai,double *ax,int &lunnz,int *luap,int *ludiag,int **luai,\
		int *rowPerm,int *rowPermInv,int *colPerm,int *colPermInv)
{
	memset(luap,0,(n+1)*sizeof(int));

	memset(ludiag,0,n*sizeof(int));

	MC21(n,ap,ai,colPerm,colPermInv);
	//MC64(n,nnz,ap,ai,ax,colPerm,colPermInv);
	MatrixPermute(n,nnz,ap,ai,ax,colPermInv);
	amd_order(n,ap,ai,rowPerm,nullptr,nullptr);
	//for(int i=0;i<n;i++)
	//{
	//	rowPerm[i]=i;
	//}
	for(int i=0;i<n;i++)
	{
		rowPermInv[rowPerm[i]]=i;
	}
	CSC2CSR(n,nnz,ap,ai,ax);
	MatrixPermute(n,nnz,ap,ai,ax,rowPerm);
	CSR2CSC(n,nnz,ap,ai,ax);
	MatrixPermute(n,nnz,ap,ai,ax,rowPerm);
	
	for(int i=0;i<n;i++)
	{
		colPerm[i]=colPermInv[rowPerm[i]];
	}
	for(int i=0;i<n;i++)
	{
		colPermInv[colPerm[i]]=i;
	}
	//符号分解
	SymbolFactor(n,ap,ai,lunnz,luap,ludiag,luai);
	
	return 1;
}

void GPLUFactorize(int n,int *ap,int *ai,double *ax,int *luap,int *ludiag,int *luai,double **luax)
{
	double *x=new double[n]{0.0};

	int start{0};
	int end{0};
	int col{0};


	int nnz=luap[n];
	if(*luax!=nullptr)
	{
		delete[] *luax;
	}
	*luax=new double[nnz]{0.0};

	for(int k=0;k<n;k++)
	{
		memset(x,0,n*sizeof(double));

		start=ap[k];
		end=ap[k+1];
		for(int i=start;i<end;i++)
		{
			x[ai[i]]=ax[i];
		}

		for(int j=luap[k];j<ludiag[k];j++)
		{
			col=luai[j];
			for(int m=ludiag[col]+1;m<luap[col+1];m++)
			{
				x[luai[m]]-=(*luax)[m]*x[col];
			}
		}

		for(int j=luap[k];j<=ludiag[k];j++)
		{
			col=luai[j];
			(*luax)[j]=x[col];
		}
		for(int j=ludiag[k]+1;j<luap[k+1];j++)
		{
			col=luai[j];
			(*luax)[j]=x[col]/x[k];
		}
	}

	if(x!=nullptr)
	{
		delete[] x;
		x=nullptr;
	}
}

void Lusolve(int n,int *luap,int *ludiag,int *luai,double *luax,double *b,int *rowPerm,int *rowPermInv,int *colPerm,int *colPermInv)
{
	int start{0};
	int end{0};

	double *x=new double[n]{0.0};
	
	for(int i=0;i<n;i++)
	{
		x[i]=b[rowPerm[i]];
	}

	//solve Ly=b
	for(int i=0;i<n-1;i++)
	{
		start=ludiag[i]+1;
		end=luap[i+1];
		for(int j=start;j<end;j++)
		{
			x[luai[j]]-=luax[j]*x[i];
		}
	}

	//solve Ux=y
	for(int i=n-1;i>=0;i--)
	{
		start=luap[i];
		end=ludiag[i];
		x[i]/=luax[end];
		for(int j=start;j<end;j++)
		{
			x[luai[j]]-=luax[j]*x[i];
		}
	}

	for(int i=0;i<n;i++)
	{
		b[i]=x[colPermInv[i]];
	}

	if(x!=nullptr)
	{
		delete[] x;
		x=nullptr;
	}
}

