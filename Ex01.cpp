
// Tian 2018-10-31

#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include<iomanip>

using namespace std;

#define N 501
// #define N 3

double vdv(double* v1, double* v2, int n);
double* mdv(double** m, double* v1, double* v2, int n);
double* vdn(double* v1, double* v2, double c, int n);
double* v2v(double* v1, double* v2, int n);
double** a2p(double a[N][N], double** p, int n);
double** aan(double** a, double c, int n);
double** ada(double** a1, double** a2, double** a3, int n);
double** storeans1(double** ans1, double* x, double lamda, int n, int c);

int powermethod1(double** a, double* x, double* lamda, int n);
int antipowermethod(double** a, double* x, double* lamda, int n);
double* GE(double** a, double* y, double* x, int n );

int main()
{
    // set the varation
    double A[N][N] = {0};
    double A2[N][N] = {0};
    int n = N;
    double** a = new double* [N];
    double** a2 = new double* [N];
    a2p(A, a, N);
    a2p(A2, a2, N);
    double X[N] = {0};
    double* x = X;
    double lamda = 0;

    // cout<<"hi";

    double ANS1[3][N+1] = {0};
    double** ans1 = new double* [N+1];
    for(int i = 0; i < 3; i++)
    {
        ans1[i] = ANS1[i];
    }

    // cout<<"hi";

    // generate the matrix
    for(int i = 0; i < N; i++ )
    {
        int k = i+1;
        a[i][i] = (1.64-0.024*(k))*sin(0.2*k)-0.64*exp(0.1/k);

        switch(k)
        {
            case 1:
            {
                a[i][1] = 0.16;
                a[i][2] = -0.064;
                break;
            }
            case 2:
            {
                a[i][2] = 0.16;
                a[i][3] = -0.064;
                a[i][0] = 0.16;
                break;
            }
            case N:
            {
                a[i][i-1] = 0.16;
                a[i][i-2] = -0.064;
                break;
            }
            case N-1:
            {
                a[i][i-1] = 0.16;
                a[i][i-2] = -0.064;
                a[i][i+1] = 0.16;
                break;
            }
            default:
            {
                a[i][i-2] = -0.064;
                a[i][i-1] = 0.16;
                a[i][i+2] = -0.064;
                a[i][i+1] = 0.16;
            }
        }   
    }    

    // read the coefficiences
    // fstream cfin("m2.txt", ios::in);
    // for(int i = 0; i < n; i++)
    // {
    //     for(int j = 0; j < n; j++)
    //     {
    //         cfin>>a[i][j];
    //     }
    // }
    // cfin.close();
    
    // solve qestion 1-1
    powermethod1(a, x, &lamda, n);

    // cout<<lamda;

    storeans1(ans1, x, lamda, n, 1);

    // solve qestion 1-2
    aan(a, -1*lamda, n);

    powermethod1(a, x, &lamda, n);

    lamda += ans1[0][N];

    // cout<<lamda;

    storeans1(ans1, x, lamda, n, 2);
        
    // solve qestion 1-3
    // ada(a, a, a2, n);   // a*a

    // for(int i = 0; i < N; i++ )
    // {
    //     for(int j = 0; j < N; j++ )
    //     {
    //         cout<<a2[i][j]<<'\t';
    //     }
    //     cout<<'\n\n';
    // }


    // powermethod1(a2, x, &lamda, n);     // the lamda max of a*a

    // cout<<lamda<<'\n';

    // ans1[2][N] = lamda;

    // aan(a2, -1*lamda, n);   // translation of a*a

    // antipowermethod(a, x, &lamda, n);     // the lamda min - lamda max of a*a

    cout<<lamda<<'\n';

    // lamda += ans1[2][N];
    // lamda = pow(lamda, 0.5);

    // cout<<lamda<<'\n';

    // lamda = ans1[2][N];

    // storeans1(ans1, x, lamda, n, 2);
        

    // print the result
    fstream cfout("result", ios::out);
    
    // cfout<<"matrix:\n";
    // for(int i = 0; i < N; i++ )
    // {
    //     for(int j = 0; j < N; j++ )
    //     {
    //         cfout<<a[i][j]<<'\t';
    //     }
    //     cfout<<'\n\n';
    // }

    cfout<<"lamda:\n";
    cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[0][N];
    cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[1][N];
    cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[2][N];
    cfout<<"\nvector:\n";
    for(int i = 0; i < n; i++)
    {
        cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[0][i];
        cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[1][i];
        cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[2][i]<<'\n';
    }    
    cfout.close();

    // release the memory
    delete a;
    delete a2;
    delete ans1;

    return 0;
}


// power method 1
int powermethod1(double** a, double* x, double* lamda, int n)
{
    // initiation
    double U[N] = {1};
    double* u = U;     
    double eta = 0;
    double Y[N] = {0};
    double* y = Y;
    double beta[2] = {0};

    double T[N] = {0};
    double* t = T;
    double b = 0;

    // cout<<"hi\n";

    for(; ; )
    {
        eta = pow(vdv(u, u, n), 0.5);
        y = vdn(u, y, 1/eta, n);
        t = mdv(a, y, t, n);
        v2v(t, u, n);
        b = beta[1];
        beta[1] = vdv(y, u, n);
        beta[0] = b;

        // cout<<beta[1]<<'\n';

        // termination judgement
        if( (abs( beta[1] - beta[0]) / abs(beta[1]) ) <= 1e-13 )
        {
            break;
        }            
    }

    v2v(u, x, n);
    (*lamda) = beta[1];
    return 0;
}


int antipowermethod(double** a, double* x, double* lamda, int n)
{
    // initiation
    double U[N] = {1};
    double* u = U;     
    double eta = 0;
    double Y[N] = {0};
    double* y = Y;
    double beta[2] = {0};

    double T[N] = {0};
    double* t = T;
    double b = 0;

    int it = 0;

    for(; ; )
    {
        eta = pow(vdv(u, u, n), 0.5);
        y = vdn(u, y, 1/eta, n);

        // t = mdv(a, y, t, n);
        // v2v(t, u, n);

        GE(a, y, u, n);     // sol equation groups to avoid reverse matrix
        it++;
        cout<<it<<'\n';

        b = beta[1];
        beta[1] = vdv(y, u, n);
        beta[0] = b;

        // cout<<beta[1]<<'\n';

        // termination judgement
        if( (abs( 1/beta[1] - 1/beta[0]) / abs(1/beta[1]) ) <= 1e-12 )
        {
            break;
        }            
    }

    v2v(u, x, n);
    (*lamda) = beta[1];
    return 0;
}


int LU(double** a, double* y, double* x, int n)
{
    double L[N][N] = {0};
    double U[N][N] = {0};

    //decompose
    for(int i = 0; i < M; i++)
    {        
        // solve U by row
        for(int j = 0; j < N; j++)
        {
            if( i > j )
            {
                U[i][j] = 0;
                continue;
            }
            else if( i == 0)
            {
                U[i][j] = a[i][j];
                continue;
            }
            else
            {
                double t = 0;                
                for(int k = 0; k < i; k++)
                {
                    t += L[i][k] * U[k][j];
                }
                U[i][j] = a[i][j] - t;
                continue;
            }
        }

        // solve L by column
        for(int j = 0; j < N; j++)
        {
            if( j == i)
            {
                L[j][i] = 1;
                continue;
            }
            else if( j < i )
            {
                L[j][i] = 0;
                continue;
            }
            else if(i == 0)
            {
                L[j][i] = a[j][i] / U[i][i];
                continue;
            }
            else
            {
                double t = 0;
                for(int k = 0; k < i; k++)
                {
                    t += L[j][k] * U[k][i];
                }
                L[j][i] = (a[j][i] - t) / U[i][i];
                continue;
            }
        }
    }

    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf("%f  ", L[i][j]);
        }
        printf("\n");
    }

    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf("%f  ", U[i][j]);
        }
        printf("\n");
    }

    // solution
    // for y
    for(int i = 0; i < M; i++)
    {
        if(i == 0)
        {
            y[i] = a[i][N];
        }
        else
        {
            double t = 0;
            for(int k = 0; k < i; k++)
            {
                t += y[k] * L[i][k];
            }
            y[i] = a[i][N] - t;
        }
    }

    // for x
    for(int i = M-1; i >= 0; i--)
    {
        if(i == (M-1))
        {
            x[i] = y[i] / U[i][i];
        }
        else
        {
            double t = 0;
            for(int k = M-1; k > i; k--)
            {
                t += x[k] * U[i][k];
            }
            x[i] = (y[i] - t) / U[i][i];
        }        
    }

    return 0;
}

double* GE(double** a, double* y, double* x, int n )
{
    // elimination
    for(int i = 0; i < n-1; i++)
    {
        for(int k = i+1; k < n; k++ )
        {
            double t = a[k][i] / a[i][i];
            
            for(int j = i; j < n; j++)
            {
                a[k][j] = a[k][j] - a[i][j] * t;
            }
            y[k] = y[k] - y[i] * t;

        }
    }

    // cout<<"elimination\n";

    // for(int i = 0; i < n; i++)
    // {
    //     for(int j = 0; j < n; j++)
    //     {
    //         printf("%f  ", a[i][j]);
    //     }
    //     printf("\n");
    // }

    // substitution
    for(int i = n-1; i >= 0; i--)
    {
        double t = 0;
        for(int j = n-1; j >= 0; j--)
        {
            if( i == (n-1) )
            {
                x[i] = y[i] / a[i][j];
                break;
            }

            else if( (j > 0) && a[i][j-1] == 0)
            {
                x[j] = ( y[i] - t ) / a[i][j];
                break;
            }
            else if( j == 0 )
            {
                x[j] = ( y[i] - t ) / a[i][j];
                break;
                
            }
            else 
            {
                t += x[j] * a[i][j];
            }
        }
    }

    // cout<<"finish\n";

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

// matrix dot vector
double* mdv(double** m, double* v1, double* v2, int n)
{
    // initiate the result
    for(int i = 0; i < n; i++)
    {
        v2[i] = 0;
    }

    // calculate
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            v2[i] += m[i][j] * v1[j];
        }
    }

    return v2;
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

// link array and pointer
double** a2p(double a[N][N], double** p, int n)
{
    for(int i = 0; i < n; i++)
    {
        // for(int j = 0; j < n; j++)
        // {
        //     p[i] = a[i];
        // }
        p[i] = a[i];
    }
    return p;
}

// matrix add
double** aan(double** a, double c, int n)
{
    for(int i = 0; i < n; i++)
    {
        a[i][i] += c;
    }
    return a;
}

double** ada(double** a1, double** a2, double** a3, int n)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                a3[i][j] += a1[i][k] * a2[j][k];
            }
        }
    }
    return a3;
}

double** storeans1(double* ans1[N+1], double* x, double lamda, int n, int c)
{
    int i = 0;
    for( ;i < n; i++)
    {
        ans1[c-1][i] = x[i];
    }
    ans1[c-1][i] = lamda;
    return ans1;
}