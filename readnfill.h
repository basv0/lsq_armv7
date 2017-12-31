/*
Here are a couple of trivial functions.
First takes data from text file having experimental points
(time/input/output) describing dynamic object, allocates memory
for the matrix of experiment X, its transpose counterpart Xt,
vector Y and fills them with the data.
Second frees memomory
*/
#ifndef _READNFILL_
#define _READNFILL_

#include <stdio.h>
#include <stdlib.h>
#include <vector>

void readNfill(FILE* yu, int M, int *pN, float ***X, float ***Xt, float **Y, float **A){
    std::vector<float> yV,uV,tV;
    char syu[80];
    float u,y,t;
    /* Takes all data into the memory in order to find out the size of matrices*/
    while(fgets(syu,80,yu)){
        sscanf(syu,"%f %f %f",&t,&u,&y);
        yV.push_back(y);
        uV.push_back(u);
        tV.push_back(t);
    }
    int N=*pN=(yV.size()-M);
    *X=new float*[N];
    *Xt=new float*[2*M];
    *A=new float[2*M];
    *Y=new float[N];
    for(int i=0;i<N;i++){
        (*X)[i]= new float[2*M];
    }
    for(int j=0;j<2*M;j++){
        (*Xt)[j]=new float[N];
    }
    for(int i=0;i<N;i++){
        (*Y)[i]=yV[i+M];
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            (*X)[i][j]=(-1.0)*yV[i+M-j-1];
            (*Xt)[j][i]=(*X)[i][j];
        }
    }
    for(int i=0; i<N;i++){
        for(int j=0;j<M;j++){
            (*X)[i][j+M]=uV[i+M-j-1];
            (*Xt)[j+M][i]=(*X)[i][j+M];
        }
    }
}

void freeMem(int M, int N, float **X, float **Xt, float *Y, float *A){
    for(int i=0;i<N;i++){
        delete[] X[i];
    }    
    for(int j=0;j<2*M;j++){
        delete[] Xt[j];
    }
    delete[] X;
    delete[] Xt;
    delete[] A;
    delete[] Y;
}

#endif
