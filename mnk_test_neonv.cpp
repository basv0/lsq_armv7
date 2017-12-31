#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>
#include "readnfill.h"
#include "mnk_neonv.h"

int M;
int N;

float **X;
float **Xt;
float *Y;
float *A;

float32x4_t **Sv;
float32x4_t *Qv;
float32x4_t **Xv;

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
    int M2 = (2*M)/Vs+(int)((2*M)%Vs!=0);
    Sv=new float32x4_t*[M2];
    Qv=new float32x4_t[N];    
    for(int i=0;i<M2;i++){
        Sv[i]=new float32x4_t[4*M];
    }
    for(int i=0;i<N;i++){
        Xv[i]=new float32x4_t[M2];
    }
    
    clock_gettime(CLOCK_REALTIME, &t1);
    mnk_neonv(N, M, M2, X, Xt, Y, A, Xv, Sv, Qv);
    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*M;i++){
        printf("%f\n",A[i]);
    }
    printf("%d %d %f\n",M,N,dt);
    freeMem(M, N, X, Xt, Y, A);
    for(int i=0;i<M2;i++){
        delete[] Sv[i];
    }
    for(int i=0;i<N;i++){
        delete[] Xv[i];
    }
    delete[] Qv;
    delete[] Sv;
    delete[] Xv;
    return 0;
}
