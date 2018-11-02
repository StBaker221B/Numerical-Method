
#include<stdio.h>

#define M 3
#define N 3

int GE(double a[M][N+1], double x[M] );

int main()
{
    // printf("hello");

    // const int M = 5, N = 5;
    double coefficients[M][N];
    double augmentation[M][N+1];
    double x[M];
    double y[M];

    // // initiation
    // for(int i = 0;i<M;i++)
    // {
    //     for(int j = 0;j<N;j++)
    //     {
    //         coefficients[i][j] = 0;
    //         augmentation[i][j] = 0;
    //     }
    //     y[i] = 0;
    //     augmentation[i][N] = 0;
    // }
    // printf("initial over \n");

    // // check
    // for(int i = 0; i < M; i++)
    // {
    //     for(int j = 0; j < N+1; j++)
    //     {
    //         printf("%f  ", augmentation[i][j]);
    //     }
    //     printf("\n");
    // }
    
    // set the coefficients
    for(int i = 0;i<M;i++)
    {
        printf("put in the coeff of line %d ", i);
        for(int j = 0;j<N;j++)
        {
            scanf("%lf", &augmentation[i][j]);
            printf("%lf \n",augmentation[i][j]);
        }
        // y[i] = 0;
        printf("put in the y of line %d ", i);
        scanf("%lf", &augmentation[i][N]);
    }
    // printf("initial over \n");
    
    // check
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            printf("%f  ", augmentation[i][j]);
        }
        printf("\n");
    }

    // solve
    GE(augmentation, x);
    for(int i = 0; i < M; i++)
    {
        printf(" %f  ", x[i]);
    }
}

int GE(double a[M][N+1], double x[M] )
{
    // elimination
    for(int i = 0; i < M-1; i++)
    {
        for(int k = i+1; k < M; k++ )
        {
            double t = a[k][i] / a[i][i];
            
            for(int j = i; j < N+1; j++)
            {
                a[k][j] = a[k][j] - a[i][j] * t;
            }
        }
    }

    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            printf("%f  ", a[i][j]);
        }
        printf("\n");
    }

    // substitution
    for(int i = M-1; i >= 0; i--)
    {
        double t = 0;
        for(int j = N-1; j >= 0; j--)
        {
            if( i == (M-1) )
            {
                x[i] = a[i][N] / a[i][j];
                break;
            }

            else if( (j > 0) && a[i][j-1] == 0)
            {
                x[j] = ( a[i][N] - t ) / a[i][j];
                break;
            }
            else if( j == 0 )
            {
                x[j] = ( a[i][N] - t ) / a[i][j];
                break;
                
            }
            else 
            {
                t += x[j] * a[i][j];
            }
        }
    }

    return 0;
}