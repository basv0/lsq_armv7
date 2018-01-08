#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <vector>
using namespace std;

#include "readnfilld.h"
#include "freemem.h"
#include "mnk_double.h"

int M;
int N;

double **X;
double **Xt;
double *Y;
double *A;
double **S;
double *Q;

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
    Q=new double[N];
    S=new double*[2*M];
    for(int j=0;j<2*M;j++){
	S[j]=new double[4*M];
    }
    
    clock_gettime(CLOCK_REALTIME, &t1);    
    
    mnk_double(N, M, X, Xt, Y, A, S, Q);
    
    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*M;i++){
        fprintf(stderr,"%f\n",A[i]);
    }
    fprintf(stderr,"t=%f\n",dt);
        printf("%d %d %f\n",M,N,dt);
    
    freeMem
    for(int j=0;j<2*M;j++){
	delete[] S[j];
    }
    delete[] Q;
    delete[] S;
    
    return 0;
}
