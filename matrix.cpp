
#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

#define N 501

int main()
{
    double a[N][N] = {0};


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

    for(int i = 0; i < N; i++ )
    {
        for(int j = 0; j < N; j++ )
        {
            cout<<a[i][j]<<'\t';
        }
        cout<<'\n';
    }

    return 0;
}