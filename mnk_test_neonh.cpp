#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "readnfill.h"
#include "mnk_neonh.h"

int M;
int N;
float **X;
float **Xt;
float **S;
float *Y;
float *Q;
float *A;

timespec t1,t2;
float dt;
FILE* yu;

int main(int argc, char *argv[]){
    if(argc<2){
        printf("least squares estimation of the dynamic model parameters\nusage: $%s n tvy\n", argv[0]);
        printf("n - model order, tvy - text file: time|input|output\n");
        return 1;
    }
    M=atoi(argv[1]);
    yu = fopen(argv[2],"r");
    if(!yu){
        printf("file %s not found\n", argv[2]);
        return 2;
    }
    readNfill(yu, M, &N, &X, &Xt, &Y, &A);
    Q=new float[N];
    S=new float*[2*M];
    for(int j=0;j<2*M;j++){
        S[j]=new float[4*M];
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    
    mnk_neonh(N, M, X, Xt, Y, A, S, Q);

    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*M;i++){
        printf("%f\n",A[i]);
    }
    printf("%d %d %f\n",M,N,dt);
    freeMem(M, N, X, Xt, Y, A);
    for(int j=0;j<2*M;j++){
        delete[] S[j];
    }
    delete[] Q;
    delete[] S;
    return 0;
}
