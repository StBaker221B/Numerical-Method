
#include<iostream>
#include<iomanip>
#include<cmath>

using namespace std;

#define N 3
#define M 3

int main()
{
    // cout<<"hello"<<endl;
    // for(int i = 0; i<5; i++)
    // {
    //     cout<<i<<'\n';
    // }
    
    double A[N][N] = {1};
    double b[2] = {0};
    double* c = b;

    cout<<A[0][1];
    cout<<A[0][0]<<endl;

    double a = 1e5 * sin(12);

    cout<<scientific<<setprecision(12)<<a<<endl;
    cout<<'\n'<<b[0];

    return 0;
}