#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<string>
#include"amd.h"

using namespace std;

int ReadMTX2CSC(string filename, int &n, int &nnz, int **ap, int **ai, double **ax);

int main()
{
	int n, nnz;
	int *ap{ nullptr };
	int	*ai{ nullptr };
	double *ax{ nullptr };

	int *perm{ nullptr };
	double control[AMD_CONTROL], info[AMD_INFO];

	if (!ReadMTX2CSC(string("test.mtx"), n, nnz, &ap, &ai, &ax))
	{
		goto END;
	}
	perm = new int[n];
	if (!perm)
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	amd_order(n, ap, ai, perm, control, info);
	amd_info(info);
	cout << endl << "perm:";
	for (int i = 0; i < n; i++)
	{
		cout << setw(10) << perm[i];
		if ((i + 1) % 10 == 0)
		{
			cout << endl;
		}
	}
	cout << endl;

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
	if (perm != nullptr)
	{
		delete[] perm;
		perm = nullptr;
	}
	
	system("pause");
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
		delete[] * ap;
	}
	*ap = new int[n + 1];
	if (!(*ap))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	if (*ai != nullptr)
	{
		delete[] * ai;
	}
	*ai = new int[nnz];
	if (!(*ai))
	{
		cout << "file \"" << __FILE__ << "\" line " << __LINE__ << ": allocate failed" << endl;
		return 0;
	}

	if (*ax != nullptr)
	{
		delete[] * ax;
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