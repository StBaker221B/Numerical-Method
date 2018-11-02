
// Tian 2018-10-31

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>

using namespace std;

#define N 501

double vdv(double* v1, double* v2, int n);
double* vdn(double* v1, double* v2, double c, int n);
double* v2v(double* v1, double* v2, int n);
double min(double a, double b);
double max(double a, double b);
double norm(double* v1, double* v2, int n);

double powermethod(double* a, double b, double c, double* x, int n);
double antipowermethod(double* a, double b, double c, double* x, double* det, int n);
double LUbounddecompose(double* a, double bb, double cc, double c[][N], int n, int r, int s);
double* LUboundsolve(double c [][N], double* B, double* x, int n, int r, int s);

int main()
{
    // set the varation

    int n = N;
    double X[N] = {0};
    double* x = X;
    double lamda1[3] = {0};
    double lamda2[N] = {0};
    double cond = 0;
    double det = 0;

    // cout<<"hi";

    double a[N] = {0};
    for(int i = 0; i < N; i++)
    {
        int k = i+1;
        a[i] = ( 1.64 - 0.024 * (k) ) * sin(0.2*k) - 0.64 * exp(0.1/k);
    }

    // cout<<"hello";

    double b = 0.16;
    double c = -0.064;
    
    // solve qestion 1-1
    lamda1[0] = powermethod(a, b, c, x, n);

    cout<<scientific<<setprecision(12)<<lamda1[0]<<'\n';

    // solve qestion 1-2

    for(int i = 0; i < N; i++)
    {
        a[i] += -1*lamda1[0];
    }

    double temp = lamda1[0];

    lamda1[1] = powermethod(a, b, c, x, n);

    lamda1[1] += temp;

    cout<<lamda1[1]<<'\n';
        
    // solve qestion 1-3

    for(int i = 0; i < N; i++)
    {
        a[i] += temp;
    }

    lamda1[2] = antipowermethod(a, b, c, x, &det, n);

    // judge lamda max and lamda min
    if( lamda1[0] > lamda1[1] )
    {
        double t = lamda1[0];
        lamda1[0] = lamda1[1];
        lamda1[1] = t;
    }
    cout<<lamda1[2]<<'\n';

    // solve question 3
    cout<<"det: "<<det<<'\n';
    cond = fabs( lamda1[0] / lamda1[2] );
    cout<<"cond2: "<<cond<<'\n';

    // solve question 2
    double s = 0;
    for(int k = 0; k < N; k++)
    {
        double u = lamda1[0] + (k+1) * ( lamda1[1] - lamda1[0] )/N;
        double a2[N] = {0};
        for(int j = 0; j < n; j++)
        {
            a2[j] = a[j] - u;
        }
        lamda2[k] = antipowermethod(a2, b, c, x, &s, n);
        lamda2[k] += u;
        cout<<(k+1)<<" :  "<<lamda2[k]<<'\n';
    }     

    // print the result
    fstream cfout("ans.txt", ios::out);
    
    cfout<<"\nWork 1\nQ1\n\n";
    cfout<<"lamda 1 : ";
    cfout<<setw(25)<<scientific<<setprecision(12)<<lamda1[0]<<'\n';
    cfout<<"lamda 501 : ";
    cfout<<setw(25)<<scientific<<setprecision(12)<<lamda1[1]<<'\n';
    cfout<<"lamda s : ";
    cfout<<setw(25)<<scientific<<setprecision(12)<<lamda1[2]<<'\n';
    
    cfout<<"\nQ2\n";
    for(int i = 0; i < N; i++)
    {
        cfout<<"lamda i "<<i+1<<" : ";
        cfout<<setw(25)<<scientific<<setprecision(12)<<lamda2[i]<<'\n';
    }    
    
    cfout<<"\nQ3\n";
    cfout<<"det(A): ";
    cfout<<setw(25)<<scientific<<setprecision(12)<<det<<'\n';
    cfout<<"cond2(A): ";
    cfout<<setw(25)<<scientific<<setprecision(12)<<cond<<'\n';

    cfout.close();

    return 0;
}

//  infinite norm
double norm(double* v1, int n)
{
    double v2 = v1[0];
    for(int i = 1; i < n; i++)
    {
        if( fabs(v1[i]) > fabs (v2) )
        {
            v2 = v1[i];
        }
    }
    return v2;
}

// power method 2
double powermethod(double* a, double b, double c, double* x, int n)
{
    // initiation
    double U[N] = {1};
    double* u = U; 

    for(int i = 0; i < n; i++)
    {
        u[i] = 1;
    }

    // double eta = 0;
    double Y[N] = {0};
    double* y = Y;
    double beta[2] = {0};

    double h0 = 0;
    double h1 = 0;

    double t = 0;

    // cout<<"hi\n";

    for(; ; )
    {
        h0 = norm(u, n);
        y = vdn(u, y, 1/fabs(h0), n);

        // calculate u
        u[0] = a[0] * y[0] + b * y[1] + c * y[2];
        u[1] = b * y[0] + a[1] * y[1] + b * y[2] + c * y[3];
        u[n-2] =  b * y[n-1] + a[n-2] * y[n-2] + b * y[n-3] + c * y[n-4];
        u[n-1] = a[n-1] * y[n-1] + b * y[n-2] + c * y[n-3];
        for(int i = 2; i < n-2; i++)
        {
            u[i] = c * ( y[i-2] + y[i+2] ) + b * ( y[i-1] + y[i+1] ) + a[i] * y[i];
        }

        h1 = norm(u, n);

        t = beta[1];
        beta[1] = ( h0/fabs(h0) )*h1;
        beta[0] = t;
        // cout<<beta[1]<<'\n';

        // termination judgement
        if( (fabs( beta[1] - beta[0]) / fabs(beta[1]) ) <= 1e-13 )
        {
            break;
        }            
    }

    v2v(u, x, n);
    return beta[1];
}

double antipowermethod(double* a, double b, double c, double* x, double* det, int n)
{
    // initiation
    double U[N] = {1};
    double* u = U;
    const int r = 2, s = 2;

    for(int i = 0; i < n; i++)
    {
        u[i] = 1;
    }

    double eta = 0;
    double Y[N] = {0};
    double* y = Y;
    double beta[2] = {0};

    double t = 0;

    double C[5][N] = {0};

    int it = 0;
    
    *det = LUbounddecompose(a, b, c, C, n, r, s);     // sol equation groups to avoid reverse matrix

    for(; ; )
    {
        eta = pow(vdv(u, u, n), 0.5);
        y = vdn(u, y, 1/eta, n);

        LUboundsolve(C, y, u, n, r, s);
        it++;
        // cout<<it<<'\n';

        t = beta[1];
        beta[1] = vdv(y, u, n);
        beta[0] = t;
        // cout<<beta[1]<<'\n';

        // termination judgement
        if( (fabs( 1/beta[1] - 1/beta[0]) / fabs(1/beta[1]) ) <= 1e-12 )
        {
            break;
        }            
    }

    v2v(u, x, n);
    return 1/beta[1];
}

// lu decompose
double LUbounddecompose(double* a, double bb, double cc, double c [][N], int n, int r, int s)
{
    // generate c
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < N; j++)
        {
            c[i][j] = 0;
        }
    }

    for(int i = 0; i < 5; i++)
    {
        switch(i)
        {
            case 0:
            {
                for(int j = 2; j < n; j++)
                {
                    c[i][j] = cc;
                }
                break;
            }

            case 1:
            {
                for(int j = 1; j < n; j++)
                {
                    c[i][j] = bb;
                }
                break;
            }

            case 2:
            {
                for(int j = 0; j < n; j++)
                {
                    c[i][j] = a[j];
                }
                break;
            }

            case 3:
            {
                for(int j = 0; j < n-1; j++)
                {
                    c[i][j] = bb;
                }
                break;
            }

            case 4:
            {
                for(int j = 0; j < n-2; j++)
                {
                    c[i][j] = cc;
                }
                break;
            }
        }
    }

    //lu decompose
    for(int k = 1; k <= n; k++)
    {
        for(int j = k; j <= min( (k+s), n ); j++)
        {
            double temp = 0;
            for(int t = max( max(1, k-r), max(k-r, j-s) ) ; t <= k-1; t++)
            {
                temp += c[k-t+s][t-1] * c[t-j+s][j-1];
            }
            c[k-j+s][j-1] = c[k-j+s][j-1] - temp;
        }

        if( k < n)
        {            
            for(int i = (k+1); i <= min((k+r), n); i++ )
            {
                double temp = 0;
                for(int t = max( max(1, (i-r)), max((i-r), (k-s) ) ); t <= k-1; t++ )
                {
                    temp += c[i-t+s][t-1] * c[t-k+s][k-1];
                }
                c[i-k+s][k-1] = (c[i-k+s][k-1] - temp) / c[s][k-1];
            }
        }
    }

    // cal the det
    double det = 1;
    for(int i = 0; i < n; i++)
    {
        det *= c[s][i];
    }

    return det;
}

// lu bond decompose
double* LUboundsolve(double c [][N], double* B, double* x, int n, int r, int s)
{
    // solve lu
    double b[N] = {0};
    for(int i = 0; i < n; i++)
    {
        b[i] = B[i];
    }

    for(int i = 2; i <= n; i++)
    {
        double temp = 0;
        for(int t = max(1, i-r); t <= i-1; t++ )
        {
            temp += c[i-t+s][t-1] * b[t-1];
        }
        b[i-1] = b[i-1] - temp;
    }

    x[n-1] = b[n-1] / c[s][n-1];
    for(int i = n-1; i >= 1; i--)
    {
        double temp = 0;
        for(int t = i+1; t <= min((i+s), n-1); t++ )
        {
            temp += c[i-t+s][t-1] * x[t-1];
        }
        x[i-1] = ( b[i-1] - temp ) / c[s][i-1];
    }

    return x;
}

// vector dot vector
double vdv(double* v1, double* v2, int n)
{
    double sum = 0;

    for(int i = 0; i < n; i++)
    {
        sum += v1[i] * v2[i];
    }

    return sum;
}

// vector dot number
double* vdn(double* v1, double* v2, double c, int n)
{
    for(int i = 0; i < n; i++)
    {
        v2[i] = v1[i] * c;
    }
    return v2;
}

// vector copy
double* v2v(double* v1, double* v2, int n)
{
    for(int i = 0; i < n; i++)
    {
        v2[i] = v1[i];
    }
    return v2;
}

double min(double a, double b)
{
    if(a > b)
        return b;
    else
        return a;
}

double max(double a, double b)
{
    if(a > b)   
        return a;
    else
        return b;
}

