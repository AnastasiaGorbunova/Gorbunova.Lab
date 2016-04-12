#include <iostream>

using namespace std;
float e = exp(1.0);

double** AllocMatrix (int );
int EnterDimention (int );
void DisplayMatrix(double** , int );
void ClearMemory (double** , int );
void MatrixA (double** , int , double );
double Sinus (double , double );
double ElementMatrixA (int , int , double );
void Matrix (double** , int );
double ElementMatrix (int , int );
void DifferenceMatrix (double** , double** , double** , int );
double MaxDifferenceMatrixElement (double** , int );

int main()
{
    const int N = 100;
    int n = EnterDimention(N);
    double** a = AllocMatrix(n);
    double eps;
    while (true)
    {
        cout << "Enter your eps: " << endl;
        cin >> eps;
        if (eps > 0 && eps <= 1) break;
        cout << "Error. Try again: " << endl;
    }
    MatrixA(a, n, eps);
    cout << "Matrix A: " << endl;
    DisplayMatrix(a, n);
    double** b = AllocMatrix(n);
    Matrix(b, n);
    cout << " Matrix B: " << endl;
    DisplayMatrix(b, n);
    double** c = AllocMatrix(n);
     DifferenceMatrix(a, b, c, n);
    cout << "Difference Matrix: " << endl;
    DisplayMatrix(c, n);
    double max = MaxDifferenceMatrixElement(c, n);
    cout << "Max element in difference matrix is  " << endl << max << endl;
    ClearMemory(a, n);
    ClearMemory(b, n);
    ClearMemory(c, n);
}

double Sinus(double x, double eps)
{
    double sum = 0;
    int k = 1;
    double n = 1;
    
    while (fabs(n) > eps)
    {
        sum += n;
        n *= (-1.0 * x * x) / ((2 * k) * (2 * k + 1));
        k += 2;
    }
    return sum;
}
void DifferenceMatrix (double** a, double** b, double** c, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[i][j] = fabs(a[i][j] - b[i][j]);
        }
    }
    return c;
}
double** AllocMatrix (int n)
{
    double** a = new double* [n];
    for (int i = 0; i < n; i++)
    {
        a[i] = new double [n];
    }
    return a;
}
void MatrixA (double** a, int n, double eps)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = ElementMyMatrix(i, j, eps);
        }
    }
    return a;
}
void Matrix (double** b, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[i][j] = ElementMatrix(i, j);
        }
    }
    return b;
}
double ElementMatrixA (int i, int j, double eps)
{
    double a = 0;
	if (i != j) a = i - j;
	else
	{
		a = (i + j + 2)*pow(e, i + j + 2) / (Sinus(2 * i + 2, eps) + 4);
	}
    return a;
}
double ElementMatrix (int i, int j)
{
    double a = 0;
	if(i != j) a = i - j;
	else
	{
		a = (i + j + 2)*pow(e, i + j + 2) / (sin(2 * i + 2) + 4);
	}
	return a;
}
double MaxDifferenceMatrixElement (double** c, int n)
{
    double max = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
			if (c[i][j] > max) max = c[i][j]; 
        }
    }
    return max;
}
int EnterDimention (int n)
{
    int k;
    while (true)
    {
        cout << "Enter size of matrix[n]*[n] 0 <= n <= " << n << ": ";
        cin >> k;
        if (k > 0 && k <= n)
            return k;
        cout << "Error!" << endl;
    }
}
void DisplayMatrix (double** a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout.width(15);
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}
void ClearMemory(double** a, int n)
{
    for (int i = 0; i < n; i++)
    {
        delete a[i];
    }
    delete []a;
}


