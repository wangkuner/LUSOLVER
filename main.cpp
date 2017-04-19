#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<string>
#include"mylu.h"

using namespace std;

int ReadMTX2CSC(string filename, int &n, int &nnz, int **ap, int **ai, double **ax);

int main()
{
	int n, nnz,lunnz;
	int *ap{ nullptr };
	int	*ai{ nullptr };
	double *ax{ nullptr };

	int *luap{nullptr};
	int *ludiag{nullptr};
	int *luai{nullptr};
	double *luax{nullptr};

	int *rowPerm{ nullptr };
	int *rowPermInv{ nullptr };
	int *colPerm{ nullptr };
	int *colPermInv{ nullptr };

	double *b{nullptr};
	double *x{nullptr};
	double *vb{nullptr};

	bool flag{true};

	if (!ReadMTX2CSC(string("test.mtx"), n, nnz, &ap, &ai, &ax))
	{
		goto END;
	}
	rowPerm = new int[4*n];
	if (!rowPerm)
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}
	rowPermInv=rowPerm+n;
	colPermInv=rowPermInv+n;
	colPerm=colPermInv+n;

	luap=new int[2*n+1];
	ludiag=luap+n+1;

	b=new double[3*n]{0.0};
	x=b+n;
	vb=x+n;

	for(int i=0;i<n;i++)
	{
		for(int j=ap[i];j<ap[i+1];j++)
		{
			b[ai[j]]+=ax[j]*(i+1.0);
		}
	}
	
	cout<<"b="<<endl;
	for(int i=0;i<n;i++)
	{
		x[i]=b[i];
		vb[i]=b[i];
		cout<<setw(10)<<x[i];
		if((i+1)%10==0)
		{
			cout<<endl;
		}
	}

	PreAnalysis(n,nnz,ap,ai,ax,lunnz,luap,ludiag,&luai,rowPerm,rowPermInv,colPerm,colPermInv);
	GPLUFactorize(n,ap,ai,ax,luap,ludiag,luai,&luax);
	Lusolve(n,luap,ludiag,luai,luax,x,rowPerm,rowPermInv,colPerm,colPermInv);

	if (!ReadMTX2CSC(string("test.mtx"), n, nnz, &ap, &ai, &ax))
	{
		goto END;
	}
	
	cout<<endl<<"colPerm="<<endl;
	for(int i=0;i<n;i++)
	{
		cout<<setw(10)<<colPerm[i];
		if((i+1)%10==0)
		{
			cout<<endl;
		}
	}

	cout<<endl<<"x="<<endl;
	for(int i=0;i<n;i++)
	{
		for(int j=ap[i];j<ap[i+1];j++)
		{
			b[ai[j]]+=ax[j]*x[i];
		}
		cout<<setw(10)<<x[i];
		if((i+1)%10==0)
		{
			cout<<endl;
		}
	}

	cout<<endl;
	for(int i=0;i<n;i++)
	{
		if(b[i]-vb[i]>1e-32&&b[i]-vb[i]<-1e-32)
		{
			cout<<i<<endl;
			flag=false;
		}
	}
	if(flag)
	{
		cout<<"right"<<endl;
	}
	else
	{
		cout<<"wrong"<<endl;
	}

END:
	if (ap != nullptr)
	{
		delete[] ap;
		ap = nullptr;
	}
	if (ai != nullptr)
	{
		delete[] ai;
		ai = nullptr;
	}
	if (ax != nullptr)
	{
		delete[] ax;
		ax = nullptr;
	}
	if (rowPerm != nullptr)
	{
		delete[] rowPerm;
		rowPerm = nullptr;
		rowPermInv=nullptr;
		colPermInv=nullptr;
		colPerm=nullptr;
	}

	if(luap!=nullptr)
	{
		delete[] luap;
		luap=nullptr;
		ludiag=nullptr;
	}
	if(luai!=nullptr)
	{
		delete[] luai;
		luai=nullptr;
	}
	if(luax!=nullptr)
	{
		delete[] luax;
		luax=nullptr;
	}

	if(b!=nullptr)
	{
		delete[] b;
		b=nullptr;
		x=nullptr;
		vb=nullptr;
	}
	
	//system("pause");
	return 0;
}

int ReadMTX2CSC(string filename, int &n, int &nnz, int **ap, int **ai, double **ax)
{
	ifstream file;
	string line, object, format, field, symmetry, row, col, value, tmp;
	istringstream head, entry;
	int i, j, ji;
	double element;
	if (filename.substr(filename.length() - 4, 4) != string(".mtx"))
	{
		cout << "file name wrong" << endl;
		return 0;
	}

	file.open(filename);
	if (!file)
	{
		cout << "open file \"" << filename << "\" failed" << endl;
		return 0;
	}

	getline(file, line);

	head.str(line);
	head >> tmp;
	head >> object;
	head >> format;
	head >> field;
	head >> symmetry;
	if (object != "matrix")
	{
		cout << "not a matrix" << endl;
		return 0;
	}
	if (format != "coordinate")
	{
		cout << "not in coo format" << endl;
		return 0;
	}

	getline(file, line);
	while (line.substr(0, 1) == "%")
	{
		getline(file, line);
	}
	entry.str(line);
	entry >> row;
	entry >> col;
	entry >> value;

	n = stoi(row);
	if (n != stoi(col))
	{
		cout << "matrix is not square" << endl;
		return 0;
	}
	nnz = stoi(value);

	if (*ap != nullptr)
	{
		delete[] *ap;
	}
	*ap = new int[n + 1];
	if (!(*ap))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	if (*ai != nullptr)
	{
		delete[] *ai;
	}
	*ai = new int[nnz];
	if (!(*ai))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	if (*ax != nullptr)
	{
		delete[] *ax;
	}
	*ax = new double[nnz];
	if (!(*ax))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	for (i = 0; i < n+1; i++)
	{
		(*ap)[i] = 0;
	}
	while (getline(file, line))
	{
		entry.str(line);
		entry.seekg(0, stringstream::beg);
		entry >> row;
		entry >> col;

		j = stoi(col);
		(*ap)[j]++;
	}

	for (i = 0; i < n; i++)
	{
		(*ap)[i + 1] += (*ap)[i];
	}

	int *count = new int[n] { 0 };
	if (!(count))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	file.clear();
	file.seekg(0, fstream::beg);
	//file.close();
	//file.open(filename);

	getline(file, line);
	while (line.substr(0, 1) == "%")
	{
		getline(file, line);
	}

	while (getline(file, line))
	{
		entry.str(line);
		entry.seekg(0, stringstream::beg);
		entry >> row;
		entry >> col;
		entry >> value;

		i = stoi(row) - 1;
		j = stoi(col) - 1;
		element = stod(value);
		
		ji = (*ap)[j] + count[j]-1;
		while (ji>=(*ap)[j]&&i < (*ai)[ji])
		{
			(*ai)[ji + 1] = (*ai)[ji];
			(*ax)[ji + 1] = (*ax)[ji];
			ji--;
		}
		(*ai)[ji+1] = i;
		(*ax)[ji+1] = element;

		count[j]++;
	}

	file.close();

	if (count != nullptr)
	{
		delete[] count;
		count = nullptr;
	}

	return 1;
}
