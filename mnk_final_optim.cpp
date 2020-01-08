/**************************************************************************
*   Copyright (C) 2019 by Basil Olonichev                                 *
*   basv0255@gmail.com                                                    *
*   v_olonichev@ksu.edu.ru                                                       *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; If not, see <http://www.gnu.org/licenses/>.  *
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>
#include <vector>

const int Vs=4; //size of vector

inline void mnk(int N, int Nr, int M, int Mr, float **X, float **Xt, float *Y, float *A, float **S, float *Q){
    
    float zero=0.0;
    float32x4_t av;
    float32x4_t bv;
    float32x4_t cv;
    float32x2_t r;
    for(int n=0;n<2*Mr;n++){
        for(int m=n;m<2*Mr;m++){
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<N/Vs;k++){
                int off=k*Vs;
                av=vld1q_f32(Xt[n]+off);
                bv=vld1q_f32(Xt[m]+off);
                cv=vmlaq_f32(cv,av,bv); 
            }
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            S[n][m]=vget_lane_f32(vpadd_f32(r,r),0);
            S[m][n]=S[n][m];
        }
        memset(S[n]+2*Mr,0,2*Mr*sizeof(float));
        S[n][n+2*Mr]=1.0;
    }
    float ratio,a,z;
    for(int i=0;i<2*Mr;i++){
        z=1.0/S[i][i];
        for(int j=0;j<2*Mr;j++){
            if(i!=j){
                ratio=S[j][i]*z;
                for(int k=0;k<4*Mr/4;k++){
                    int off=k*Vs;
                    av=vld1q_f32(S[i]+off);
                    bv=vld1q_f32(S[j]+off);
                    bv=vmlsq_n_f32(bv,av,ratio);
                    vst1q_f32(S[j]+off,bv);
                }       
            }
        }
    }
    for(int i=0;i<2*Mr;i++){
        a=1.0/S[i][i];
        for(int j=0;j<Mr;j++){
            int off=j*Vs;
            av=vld1q_f32(S[i]+off);
            av=vmulq_n_f32(av,a);
            vst1q_f32(S[i]+off,av);     
        }
    }
    for(int n=0;n<2*Mr;n++){
        for(int m=0;m<Nr;m++){
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<2*M/Vs;k++){
                int off=k*Vs;
                av=vld1q_f32(S[n]+2*Mr+off);
                bv=vld1q_f32(X[m]+off);
                cv=vmlaq_f32(cv,av,bv);         
            }
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            Q[m]=vget_lane_f32(vpadd_f32(r,r),0);
        }
        cv=vld1q_dup_f32(&zero);
        for(int m=0;m<N/Vs;m++){
            int off=m*Vs;
            av=vld1q_f32(Q+off);
            bv=vld1q_f32(Y+off);
            cv=vmlaq_f32(cv,av,bv);         
        }
        r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
        A[n]=vget_lane_f32(vpadd_f32(r,r),0);       
    }
}

int M,Mr;    //nearest to multiple of 4 and real size
int N,Nr;
float **X;
float **Xt;
float **S;
float *Y;
float *Q;
float *A;

timespec t1,t2;
float dt;

int main(int argc, char *argv[]){
    if(argc<3){	
        printf("usage $%s f M\n",argv[0]);
        printf("f-text file with the experimant points: time input output, M-objects order\n");
        return 0;
    }
    FILE* yu = fopen(argv[1],"r");
    if(!yu){
        printf("file %s not found\n", argv[1]);
        return 1;
    }
    Mr=atoi(argv[2]);    
    std::vector<float> yV,uV,tV;
    char syu[80];
    float u,y,t;
    while(fgets(syu,80,yu)){
        sscanf(syu,"%f %f %f",&t,&u,&y);
        yV.push_back(y);
        uV.push_back(u);
        tV.push_back(t);
    }
    fclose(yu);
    Nr=(yV.size()-Mr);
    
    int D=Nr/Vs;
    int R=Nr%Vs;
    D+=(R!=0);
    int N=D*Vs;
    int Ds=(2*Mr)/Vs;
    int Rs=(2*Mr)%Vs;
    Ds+=(Rs!=0);
    int M=Ds*Vs/2;
    
    X=new float*[N];
    Xt=new float*[2*M];
    Q=new float[N];
    S=new float*[2*Mr];
    A=new float[2*Mr];
    Y=new float[N];
    
    for(int i=0;i<N;i++){
        X[i]=new float[2*M];
    }
    
    for(int i=0;i<2*M;i++){
        Xt[i]=new float[N];
    }
    for(int i=0;i<2*Mr;i++){
        S[i]=new float[2*Mr+2*M];
    }
    
    for(int i=0;i<Nr;i++){
        Y[i]=yV[i+Mr];
    }
    for(int i=Nr;i<N;i++){ 
        Y[i]=0.0;//left-over
        Q[i]=0.0;
    }
    for(int i=0;i<Nr;i++){
        for(int j=0;j<Mr;j++){
            X[i][j]=(-1.0)*yV[i+Mr-j-1];
            Xt[j][i]=X[i][j];
        }
    }
    for(int i=0; i<Nr;i++){
        for(int j=0;j<Mr;j++){
            X[i][j+Mr]=uV[i+Mr-j-1];
            Xt[j+Mr][i]=X[i][j+Mr];
        }
        for(int j=2*Mr;j<2*M;j++){
            X[i][j]=0.0;//left-over
            Xt[j][i]=X[i][j];
        }
    }
    for(int i=Nr;i<N;i++){
        for(int j=0;j<2*M;j++){
            X[i][j]=0.0;//left-over
            Xt[j][i]=X[i][j];
        }
    }
    
    clock_gettime(CLOCK_REALTIME, &t1);
    
    mnk(N, Nr, M, Mr, X, Xt, Y, A, S, Q);
    
    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*Mr;i++){
        fprintf(stderr,"%f\n",A[i]);
    }
    fprintf(stderr,"t=%f\n",dt);
    printf("%d %d %f\n",Mr,Nr,dt);

    for(int i=0;i<N;i++){
        delete[] X[i];
    }
    
    for(int j=0;j<2*M;j++){        
        delete[] Xt[j];
    }
    for(int j=0;j<2*Mr;j++){        
        delete[] S[j];
    }
    delete[] X;
    delete[] Xt;
    delete[] Q;
    delete[] S;
    delete[] A;
    delete[] Y;
    
    return 0;
}
