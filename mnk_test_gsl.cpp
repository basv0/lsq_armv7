#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include </usr/local/include/gsl/gsl_matrix.h>
#include </usr/local/include/gsl/gsl_linalg.h>
#include </usr/local/include/gsl/gsl_permutation.h>
#include <vector>
//g++ -O2 -o mnk_gsl ./mnk_gsl.cpp -lgsl -lgslcblas -lrt
FILE* yu;
int M=3;
int d=0;
int N;

std::vector<double> yV,uV,tV;
char syu[80];
gsl_matrix *F;
gsl_vector *T,*R,*S,*Y;
timespec t1,t2;
float dt;

int main(int argc, char* argv[]){
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
    std::vector<double> yV,uV,tV;                    
    char syu[80];                                    
    double u,y,t;                                    
    while(fgets(syu,80,yu)){                         
	sscanf(syu,"%lf %lf %lf",&t,&u,&y);          
	yV.push_back(y);                             
	uV.push_back(u);                             
	tV.push_back(t);                             
    }                                                
    N=(yV.size()-M);                                 
    
    F=gsl_matrix_calloc(N,2*M);
    T = gsl_vector_calloc(2*M);
    R = gsl_vector_calloc(2*M);
    S = gsl_vector_calloc(N);
    Y = gsl_vector_calloc(N);

    for(int i=0; i<N; i++){
	gsl_vector_set(Y,i,yV[i+M]);
    }

    for(int i=0; i<N; i++){
        for(int j=0;j<M;j++){
	    gsl_matrix_set(F,i,j,(-1.0)*yV[i+M-j-1]);
        }
    }

    for(int i=0; i<N;i++){
        for(int j=0;j<M;j++){
	    gsl_matrix_set(F,i,j+M,uV[i+M-j-1]);
        }
    }
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);    
    gsl_linalg_QR_decomp(F,T);
    gsl_linalg_QR_lssolve(F,T,Y,R,S);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int j=0; j<M+M; j++){
        fprintf(stderr,"%f\n",gsl_vector_get(R,j));
    }
    fprintf(stderr,"t=%f\n",dt);
    printf("%d %d %f\n",M,N,dt);    
    gsl_matrix_free(F);
    gsl_vector_free(T);
    gsl_vector_free(R);
    gsl_vector_free(S);
    gsl_vector_free(Y);
    
    return 0;
}
