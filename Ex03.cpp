
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
double min(double a, double b);
double max(double a, double b);
double norm(double* v1, double* v2, int n);


double** storeans1(double** ans1, double* x, double lamda, int n, int c);

double powermethod(double* a, double b, double c, double* x, int n);
double antipowermethod(double* a, double b, double c, double* x, int n);
double* LU(double* a, double b, double c, double* B, double* x, int n);
double* LUbound(double* a, double b, double c, double* B, double* x, int n);

int main()
{
    // set the varation

    // double A[N][N] = {0};
    // double A2[N][N] = {0};
    int n = N;
    // double** a = new double* [N];
    // double** a2 = new double* [N];
    // a2p(A, a, N);
    // a2p(A2, a2, N);
    double X[N] = {0};
    double* x = X;
    double lamda = 0;

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

    // double ANS1[3][N+1] = {0};
    // double** ans1 = new double* [N+1];
    // for(int i = 0; i < 3; i++)
    // {
    //     ans1[i] = ANS1[i];
    // }

    // cout<<"hi";
    
    // solve qestion 1-1
    lamda = powermethod(a, b, c, x, n);

    cout<<scientific<<setprecision(12)<<lamda<<'\n';
    
    // for(int i = 0; i < N; i++)
    // {
    //     a[i] += -10;
    // }

    // lamda = antipowermethod(a, b, c, x, n);

    // cout<<lamda<<'\n';



    // storeans1(ans1, x, lamda, n, 1);

    // solve qestion 1-2
    // aan(a, -1*lamda, n);

    for(int i = 0; i < N; i++)
    {
        a[i] += -1*lamda;
    }

    double temp = lamda;

    lamda = powermethod(a, b, c, x, n);

    lamda += temp;

    cout<<lamda<<'\n';

    // powermethod1(a, x, &lamda, n);

    // lamda += ans1[0][N];

    // cout<<lamda;

    // storeans1(ans1, x, lamda, n, 2);
        
    // solve qestion 1-3

    for(int i = 0; i < N; i++)
    {
        a[i] += temp;
    }

    lamda = antipowermethod(a, b, c, x, n);

    cout<<lamda;

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

    // cout<<lamda<<'\n';

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

    // cfout<<"lamda:\n";
    // cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[0][N];
    // cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[1][N];
    // cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[2][N];
    // cfout<<"\nvector:\n";
    // for(int i = 0; i < n; i++)
    // {
    //     cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[0][i];
    //     cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[1][i];
    //     cfout<<setw(25)<<scientific<<setprecision(12)<<ans1[2][i]<<'\n';
    // }    
    // cfout.close();

    // // release the memory
    // delete a;
    // delete a2;
    // delete ans1;

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

    // double t[N] = {0};
    // double* t = T;
    double t = 0;

    // cout<<"hi\n";

    for(; ; )
    {
        // eta = pow(vdv(u, u, n), 0.5);
        // y = vdn(u, y, 1/eta, n);
        // t = mdv(a, y, t, n);
        // v2v(t, u, n);

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
        // beta[1] = vdv(y, u, n);
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
    // (*lamda) = beta[1];
    return beta[1];
}


double antipowermethod(double* a, double b, double c, double* x, int n)
{
    // initiation
    double U[N] = {1};
    double* u = U;

    for(int i = 0; i < n; i++)
    {
        u[i] = 1;
    }


    double eta = 0;
    double Y[N] = {0};
    double* y = Y;
    double beta[2] = {0};

    // double T[N] = {0};
    // double* t = T;
    double t = 0;

    int it = 0;

    for(; ; )
    {
        eta = pow(vdv(u, u, n), 0.5);
        y = vdn(u, y, 1/eta, n);

        // t = mdv(a, y, t, n);
        // v2v(t, u, n);

        u = LUbound(a, b, c, y, u, n);     // sol equation groups to avoid reverse matrix
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
    // (*lamda) = beta[1];
    return 1/beta[1];
}

// lu decompose
double* LU(double* a, double b, double c, double* B, double* x, int n)
{
    double l[N][N] = {1};
    double u[N][N] = {0};
    double y[N] = {0};

    // for u and l
    for(int k = 1; k <= n; k++ )
    {
        // u
        for(int j = k; j <= min((k+2), n); j++)
        {
            if(j == 1)
            {
                u[0][0] = a[0];
                continue;
            }
            else
            {
                double sum = 0;

                for(int t = max( max(1,(k-2)), max((k-2),(j-2)) ); t <= k-1; t++)
                {
                    // cout<<t;
                    sum += l[k-1][t-1] * u[t-1][j-1];
                }
                // u[k-1][j-1] = a[k-1][j-1] - sum;
                if(j == k)
                {
                    u[k-1][j-1] = a[k-1] - sum;
                    continue;
                }
                else if( j == k+1)
                {
                    u[k-1][j-1] = b - sum;
                    continue;
                }
                else if( j == k+2)
                {
                    u[k-1][j-1] = c - sum;
                    continue;
                }
            }
        }
    
        // l
        for(int i = k+1; i <= min((k+2), n); i++)
        {
            // if(i == 1)
            // {
            //     l[0][0] = 1;
            //     continue;
            // }

            if(k == n)
            {
                l[n-1][n-1] = 1;
                continue;
            }
            else 
            {
                double sum = 0;
                for(int t = max( max(1,(i-2)), max((k-2),(i-2)) ); t <= k-1; t++)
                {
                    sum += l[i-1][t-1] * u[t-1][k-1];
                }
                // l[i-1][k-1] = ( a[i-1][k-1] - sum ) / u[k-1][k-1];
                if( i == k+1)
                {
                    l[i-1][k-1] = ( b - sum ) / u[k-1][k-1];
                    continue;
                }
                else if( i == k+2 )
                {
                    l[i-1][k-1] = ( c - sum ) / u[k-1][k-1];
                    continue;
                }
            }
        }
        
    }
    
    // solve y
    y[0] = B[0];
    double temp = 0;
    for(int i = 2; i <= n; i++)
    {
        for(int t = max(1,i-2); t <= i-1; t++ )
        {
            temp += l[i-1][t-1] * y[t-1];
        }
        y[i-1] = B[i-1] - temp;
    }

    //solve x
    temp = 0;
    x[n-1] = y[n-1] / u[n-1][n-1];
    for(int i = n-1; i >= 1; i--)
    {
        for(int t = i+1; t <= min(i+2, n); t++ )
        {
            temp += u[i-1][t-1] * x[t-1];
        }
        x[i-1] = ( y[i-1] - temp ) / u[i-1][i-1] ;
    }

    return x;
}

// lu bond decompose
double* LUbound(double* a, double bb, double cc, double* B, double* x, int n)
{
    // generate c
    double c[5][N] = {0};
    // for(int i = 0; i < 5; i++)
    // {
    //     for(int j = 0; j < N; j++)
    //     {
    //         c[i][j] = 0;
    //     }
    // }

    // for(int i = 0; i < 5; i++)
    // {
    //     for(int j = 0; j < N; j++)
    //     {
    //         cout<<c[i][j];
    //     }
    //     cout<<'\n';
    // }

    const int s = 2, r = 2;
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
        // if( (i == 4) || (i == 0) )
        // {
        //     for(int j = s-i)
        // }
        // for(int j =  ; j < n; j++)
        // {
        //     if( (i-s) < 0 )
        //     {
        //         c[i][]
        //     }
        // }
        /*
        for(int j = 0; j < n; j++)
        {
            if( j == i)
            {
                c[s][j] = a[i];
                break;
            }
            else if( j == i-1)
            {
                c[1+s][j] = bb;
                break;
            }
            else if( j == i-2)
            {
                c[2+s][j] = cc;
                break;
            }
            else if( j == i+1)
            {
                c[s-1][j] = bb;
                break;
            }
            else if( j == i+2)
            {
                c[s-2][j] = cc;
                break;
            }            
        }
        */
    }

    // for(int i = 0; i < N; i++)
    // {
    //     for(int j = 0; j < 5; j++)
    //     {
    //         cout<<c[j][i];
    //     }
    //     cout<<'\n';
    // }

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
                // cout<<k;
                c[i-k+s][k-1] = (c[i-k+s][k-1] - temp) / c[s][k-1];
                // cout<<i-k+s+1<<'\n';
            }
            // cout<<i<<'\n';
            // cout<<k<<'\n';
        }
    }

    // show c
    // for(int i = 0; i < 5; i++)
    // {
    //     for(int j = 0; j < N; j++)
    //     {
    //         cout<<c[i][j];
    //     }
    //     cout<<'\n';
    // }

    // solve lu
    double b[N] = {0};
    // b = v2v(B, b, n);
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

double min(double a, double b)
{
    if(a > b)
    {
        return b;
    }
    else
    {
        return a;
    }
}

double max(double a, double b)
{
    if(a > b)
    {
        return a;
    }
    else
    {
        return b;
    }
}