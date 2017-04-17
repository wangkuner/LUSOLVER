#include<iostream>
#include<cstring>
#include<iterator>

#include"mylu.h"
#include"mc64.h"
#include"amd.h"

using std::cout;
using std::endl;
using std::vector;



void FormatMatrixIntoCSR(Matrix2D<double> matA,long &n,long &nnz,long **ap,long **ai,double **ax)
{
	unsigned long num{0};

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
	*ap=new long[n+1]{0};

	for(unsigned long i=0;i<matA.row;i++)
	{
		for(unsigned long j=0;j<matA.col;j++)
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
	*ai=new long[nnz]{0};

	num=0;

	for(unsigned long i=0;i<matA.row;i++)
	{
		for(unsigned long j=0;j<matA.col;j++)
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

void CSR2CSC(long n,long nnz,long *ap,long *ai,double *ax)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配CSC存储空间"<<endl;
		exit(0);
	}

	size_t size=0;

	int *nonzeropercol=new int[n]{0};

	long *cap=new long[n+1]{0};
	long *cai=new long[nnz]{0};
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

	size=(n+1)*sizeof(long);
	memcpy(ap,cap,size);
	size=nnz*sizeof(long);
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

void CSC2CSR(long n,long nnz,long *ap,long *ai,double *ax)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配CSR存储空间"<<endl;
		exit(0);
	}

	size_t size=0;

	int *nonzeropercol=new int[n]{0};

	long *cap=new long[n+1]{0};
	long *cai=new long[nnz]{0};
	double *cax=new double[nnz]{0};

	size=(n+1)*sizeof(long);
	memcpy(cap,ap,size);
	memset(ap,0,size);
	size=nnz*sizeof(long);
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

void MatrixPermute(long n,long nnz,long *ap,long *ai,double *ax,long *perm)
{
	if(ap==nullptr||ai==nullptr||ax==nullptr)
	{
		cout<<"未分配存储空间"<<endl;
		exit(0);
	}

	size_t size=0;
	long p=0;

	long *tap=new long[n+1]{0};
	long *tai=new long[nnz]{0};
	double *tax=new double[nnz]{0};

	size=(n+1)*sizeof(long);
	memcpy(tap,ap,size);
	size=(nnz)*sizeof(long);
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

int DFS(long n,long col,long *lucap,long *ludiag,std::vector<long> lucai,char *line)
{
	/*递归方法速度太慢
	if(line[colki]==VISITED)
	{
		return 0;
	}

	long cstart{0};
	long cend{0};

	line[colki]=VISITED;

	if(lucap[colki+1]-ludiag[colki]>1)
	{
		cstart=ludiag[colki]+1;
		cend=lucap[colki+1];
		for(long i=cstart;i<cend;i++)
		{
			DFS(lucai[i],lucap,ludiag,lucai,line);
		}
	}
	else
	{
		return 0;
	}
	*/
	long coli{ludiag[col]+1};
	vector<long> stack(2*n);
	int top=-1;

	
	line[col]=VISITED;

	while(coli<lucap[col+1] || (!(top==-1)))
	{
		if(coli<lucap[col+1])
		{
			while(coli<lucap[col+1])
			{
				if(line[lucai[coli]]!=VISITED)
				{
					stack[++top]=col;
					stack[++top]=coli;
					col=lucai[coli];
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
			while(coli<lucap[col+1])
			{
				if(line[lucai[coli]]!=VISITED)
				{
					//stack.push_back(col);
					//stack.push_back(coli);
					stack[++top]=col;
					stack[++top]=coli;
					col=lucai[coli];
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

int SymbolPrediction(long n,long k,long *cap,long *cai,long *lucap,long *ludiag,std::vector<long> lucai,char *nonzero)
{
	char *line=new char[n]{0};

	long cstart=cap[k];
	long cend=cap[k+1];
	long ln=0;

	for(int i=cstart;i<cend;i++)
	{
		DFS(n,cai[i],lucap,ludiag,lucai,line);
		
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

int SymbolFactor(long n,long *ap,long *ai,long &lunnz,long *luap,long *ludiag,long **luai)
{
	vector<long>tmpluai;

	char *nonzero=new char[n]{0};

	if(!nonzero)
	{
		cout<<"file \""<<__FILE__<<"\" line "<<__LINE__<<":"<<"allocating memory failed"<<endl;
		return 0;
	}

	//symbol prediction, get the parttern of LU
	for(int k=0;k<n;k++)
	{
		SymbolPrediction(n,k,cap,*cai,lucap,ludiag,lucai,nonzero);

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
		luap[k+1]=lucap[k]+len;
		ludiag[k]=lucap[k]+diag-1;
	}

	lunnz=tmpluai.size();

	if(*luai!=nullptr)
	{
		delete[] *luai;
	}
	*luai=new long[lunnz];
	if(!(*luai))
	{
		cout<<"file \""<<__FILE__<<"\" line "<<__LINE__<<":"<<"allocating memory failed"<<endl;
		return 0;
	}

	for(long i=0;i<lunnz;i++)
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

int PreAnalysis(long n,long nnz,long *ap,long *ai,double *ax,long &lunnz,long *luap,long *ludiag,long **luai,\
		long *rowPerm,long *rowPermInv,long *colPerm,long *colPermInv)
{
	if((end(luap)-begin(luap))!=n+1)
	{
		cout<<"memory space of luap is insuffient"<<endl;
		return 0;
	}
	memset(luap,0,(n+1)*sizeof(long));

	if((end(ludiag)-begin(ludiag))!=n)
	{
		cout<<"memory space of ludiag is insuffient"<<endl;
		return 0;
	}
	memset(ludiag,0,n*sizeof(long));

	if((end(rowPerm)-begin(rowPerm))!=n)
	{
		cout<<"memory space of rowPerm is insuffient"<<endl;
		return 0;
	}

	if((end(rowPermInv)-begin(rowPermInv))!=n)
	{
		cout<<"memory space of rowPermInv is insuffient"<<endl;
		return 0;
	}

	if((end(colPerm)-begin(colPerm))!=n)
	{
		cout<<"memory space of colPerm is insuffient"<<endl;
		return 0;
	}

	if((end(colPermInv)-begin(colPermInv))!=n)
	{
		cout<<"memory space of colPermInv is insuffient"<<endl;
		return 0;
	}

	MC64(n,nnz,ap,ai,ax,colPerm,colPermInv);
	MatrixPermute(n,nnz,ap,ai,ax,colPermInv);
	amd_l_order(n,ap,ai,rowPerm,nullptr,nullptr);
	for(long i=0;i<n;i++)
	{
		rowPermInv[rowPerm[i]]=i;
	}
	CSC2CSR(n,nnz,ap,ai,ax);
	MatrixPermute(n,nnz,ap,ai,ax,rowPermInv);
	CSR2CSC(n,nnz,ap,ai,ax);
	//符号分解
	SymbolFactor(n,ap,ai,lunnz,luap,ludiag,luai);

	return 1;
}

void GPLUFactorize(long n,long &nnz,long *cap,long *cdiag,long **cai,double **cax)
{
	vector<long> lucap(n+1,0);
	vector<long> lucai;
	vector<long> ludiag(n,0);
	
	double *lucax{nullptr};
	long *tai{nullptr};

	vector<char> nonzero(n,0);

	size_t size=n*sizeof(double);
	double *x=new double[n]{0.0};

	long cstart{0};
	long cend{0};
	long col{0};

	int len{0};
	int diag{0};

	//symbol prediction, get the parttern of LU
	for(int k=0;k<n;k++)
	{
		for(auto &t:nonzero)
		{
			t=0;
		}
		SymbolPrediction(n,k,cap,*cai,lucap,ludiag,lucai,nonzero);
		/*
		for(auto &t:x)
		{
			t=0.0;
		}
		cstart=cap[k];
		cend=cap[k+1];
		for(long i=cstart;i<cend;i++)
		{
			x[(*cai)[i]]=(*cax)[i];
		}

		for(int j=0;j<k;j++)
		{
			if(nonzero[j]==VISITED)
			{
				for(int m=ludiag[j]+1;m<lucap[j+1];m++)
				{
					x[lucai[m]]-=lucax[m]*x[j];
				}
			}
		}
		*/

		len=0;
		diag=0;
		for(int i=0;i<n;i++)
		{
			if(i<=k)
			{
				//if(x[i]<-1e-32||x[i]>1e-32)
				if(nonzero[i]!=0)
				{
					len++;
					diag++;
					lucai.push_back(i);
					//lucax.push_back(x[i]);
				}
			}
			else
			{
				//x[i]/=x[k];
				//if(x[i]<-1e-32||x[i]>1e-32)
				if(nonzero[i]!=0)
				{
					len++;
					lucai.push_back(i);
					//lucax.push_back(x[i]);
				}
			}
		}
		lucap[k+1]=lucap[k]+len;
		ludiag[k]=lucap[k]+diag-1;
	}
	//end

	nnz=lucap[n];
	lucax=new double[nnz]{0.0};
	tai=new long[nnz]{0};

	for(int k=0;k<n;k++)
	{
		memset(x,0,size);

		cstart=cap[k];
		cend=cap[k+1];
		for(long i=cstart;i<cend;i++)
		{
			x[(*cai)[i]]=(*cax)[i];
		}

		for(int j=lucap[k];j<ludiag[k];j++)
		{
			col=lucai[j];
			for(int m=ludiag[col]+1;m<lucap[col+1];m++)
			{
				x[lucai[m]]-=lucax[m]*x[col];
			}
		}

		for(int j=lucap[k];j<=ludiag[k];j++)
		{
			col=lucai[j];
			tai[j]=col;
			lucax[j]=x[col];
		}
		for(int j=ludiag[k]+1;j<lucap[k+1];j++)
		{
			col=lucai[j];
			tai[j]=col;
			lucax[j]=x[col]/x[k];
		}
	}

	/*
	if(lucap[n]!=lucai.size())
	{
		cout<<"LU factor wrong"<<endl;
	}

	if(nnz!=lucap[n])
	{
		nnz=lucap[n];
		delete[] *cai;
		delete[] *cax;

		*cai=new long[nnz]{0};
		*cax=new double[nnz]{0.0};

	}*/

	cap[0]=lucap[0];
	for(int i=0;i<n;i++)
	{
		cap[i+1]=lucap[i+1];
		cdiag[i]=ludiag[i];
	}

	if(*cai!=nullptr)
	{
		delete[] *cai;
	}
	if(*cax!=nullptr)
	{
		delete[] *cax;
	}
	*cai=tai;
	*cax=lucax;
	/*
	for(int i=0;i<nnz;i++)
	{
		(*cai)[i]=lucai[i];
		(*cax)[i]=lucax[i];
	}*/
}
void MyNicsluSolve(long n,long *cap,long *cdiag,long *cai,double *cax,double *b,long *rowPerm,long *rowPermInv,long *colPerm,long *colPermInv)
{
	long cstart{0};
	long cend{0};

	double *x=new double[n]{0.0};
	
	for(int i=0;i<n;i++)
	{
		x[i]=nicslu->row_scale[nicslu->row_perm[i]]*b[nicslu->row_perm[i]];
	}

	//solve Ly=b
	for(int i=0;i<n-1;i++)
	{
		cstart=cdiag[i]+1;
		cend=cap[i+1];
		for(long j=cstart;j<cend;j++)
		{
			x[cai[j]]-=cax[j]*x[i];
		}
	}

	//solve Ux=y
	for(int i=n-1;i>=0;i--)
	{
		cstart=cap[i];
		cend=cdiag[i];
		x[i]/=cax[cend];
		for(long j=cstart;j<cend;j++)
		{
			x[cai[j]]-=cax[j]*x[i];
		}
	}

	for(int i=0;i<n;i++)
	{
		b[i]=nicslu->col_scale_perm[nicslu->col_perm_inv[i]]*x[nicslu->col_perm_inv[i]];
	}

	if(x!=nullptr)
	{
		delete[] x;
		x=nullptr;
	}
}
//__global__ void LUFactorKernel(long n,long *ap,long *ai,double *ax,long* luap,long *ludiag,long *luai,double *luax)
//{
//}
//
//__global__ void SolveKernel(long n,long* luap,long *ludiag,long *luai,double *luax,double *b)
//{
//}

