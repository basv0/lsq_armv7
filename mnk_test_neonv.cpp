#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>

#include <vector>
using namespace std;

#include "readnfill.h"
#include "freemem.h"
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
FILE *yu;

int main(int argc, char *argv[]){
    if(argc<3){
        printf("least squares estimation of the dynamic model parameters\nusage: $%s n tvy\n", argv[0]);
        printf("n - model order, tvy - text file: time|input|output\n");
        return 1;
    }

    readNfill
    fclose(yu);
    
    int M2 = (2*M)/Vs+(int)((2*M)%Vs!=0);
    delete[] A;
    A=new float[M2*Vs];   
    Sv=new float32x4_t*[M2];
    Qv=new float32x4_t[N];
    Xv=new float32x4_t*[N];
    
    for(int i=0;i<M2;i++){
        Sv[i]=new float32x4_t[4*M];
    }
    
    
    for(int i=0;i<N;i++){
        Xv[i]=new float32x4_t[M2];
    }
    
    
    int R = (2*M)%Vs;
    for(int i=0;i<N;i++){
        for(int j=0;j<M2-1;j++){
            Xv[i][j]=vld1q_f32(X[i]+j*Vs);
        }
        if(!R){
            Xv[i][M2-1]=vld1q_f32(X[i]+(M2-1)*Vs);
        }
        else{
            float tail[Vs];
            memset(tail,0,Vs*sizeof(float));
            memcpy(tail,X[i]+(M2-1)*Vs,R*sizeof(float)); 
            Xv[i][M2-1]=vld1q_f32(tail);
        }
    }
    clock_gettime(CLOCK_REALTIME, &t1);    
    
    mnk_neonv(N, M, M2, X, Xt, Y, A, Xv, Sv, Qv);
    
    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*M;i++){
        fprintf(stderr,"%f\n",A[i]);
    }
    fprintf(stderr,"t=%f\n",dt);
        printf("%d %d %f\n",M,N,dt);
    
    freeMem
    
    for(int i=0;i<N;i++){
        delete[] Xv[i];
    }
    
    for(int i=0;i<M2;i++){
        delete[] Sv[i];
    }
    
    delete[] Xv;
    delete[] Qv;
    delete[] Sv;
    
    return 0;
}
