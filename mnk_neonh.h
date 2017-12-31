/**************************************************************************
*   Copyright (C) 2010 by Basil Olonichev                                 *
*   basv0255@gmail.com                                                    *
*   ovv@kstu.edu.ru                                                       *
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

#ifndef _MNK_NEONH_
#define _MNK_NEONH_

#include <arm_neon.h>
#include <string.h>

const int Vs=128/sizeof(float)/8;

inline void mnk_neonh(int N, int M, float **X, float **Xt, float *Y, float *A, float **S, float *Q){
    int D=N/Vs;
    int R=N%Vs;
    float zero=0.0;
    float32x4_t av;
    float32x4_t bv;
    float32x4_t cv;
    float32x2_t r;
    float tail_a[Vs];
    float tail_b[Vs];    
    memset(tail_a,0,Vs*sizeof(float));
    memset(tail_b,0,Vs*sizeof(float));
    for(int n=0;n<2*M;n++){
        memcpy(tail_a,Xt[n]+D*Vs,R*sizeof(float)); 
        for(int m=n;m<2*M;m++){
            memcpy(tail_b,Xt[m]+D*Vs,R*sizeof(float)); 
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<D;k++){
                int off=k*Vs;
                av=vld1q_f32(Xt[n]+off);
                bv=vld1q_f32(Xt[m]+off);
                cv=vmlaq_f32(cv,av,bv);         
            }
            av=vld1q_f32(tail_a);
            bv=vld1q_f32(tail_b);
            cv=vmlaq_f32(cv,av,bv); 
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            S[n][m]=vget_lane_f32(vpadd_f32(r,r),0);
            if(m!=n) S[m][n]=S[n][m];
        }
        memset(S[n]+2*M,0,2*M*sizeof(float));
        S[n][n+2*M]=1.0;
    }
    float ratio,a,z;
    for(int i=0;i<2*M;i++){
        z=1.0/S[i][i];
        for(int j=0;j<2*M;j++){
            if(i!=j){
                ratio=S[j][i]*z;
                for(int k=0;k<M;k++){
                    int off=k*Vs;
                    av=vld1q_f32(S[i]+off);
                    bv=vld1q_f32(S[j]+off);
                    bv=vmlsq_n_f32(bv,av,ratio);
                    vst1q_f32(S[j]+off,bv);
                }       
            }
        }
    }
    for(int i=0;i<2*M;i++){
        a=1.0/S[i][i];
        for(int j=0;j<M;j++){
            int off=j*Vs;
            av=vld1q_f32(S[i]+off);
            av=vmulq_n_f32(av,a);
            vst1q_f32(S[i]+off,av);     
        }
    }
    D = (2*M)/Vs;
    R = (2*M)%Vs;
    memset(tail_a,0,Vs*sizeof(float));
    memset(tail_b,0,Vs*sizeof(float));
    for(int n=0;n<2*M;n++){
        memcpy(tail_a,S[n]+2*M+D*Vs,R*sizeof(float));   
        for(int m=0;m<N;m++){
            memcpy(tail_b,X[m]+D*Vs,R*sizeof(float)); 
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<D;k++){
                int off=k*Vs;
                av=vld1q_f32(S[n]+2*M+off);
                bv=vld1q_f32(X[m]+off);
                cv=vmlaq_f32(cv,av,bv);         
            }
            av=vld1q_f32(tail_a);
            bv=vld1q_f32(tail_b);
            cv=vmlaq_f32(cv,av,bv); 
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            Q[m]=vget_lane_f32(vpadd_f32(r,r),0);       
        }
        int D=N/Vs;
        int R=N%Vs;
        float tail_c[Vs];
        float tail_d[Vs];
        cv=vld1q_dup_f32(&zero);
        memcpy(tail_c,Q+D*Vs,R*sizeof(float));    
        memcpy(tail_d,Y+D*Vs,R*sizeof(float));    
        for(int m=0;m<D;m++){
            int off=m*Vs;
            av=vld1q_f32(Q+off);
            bv=vld1q_f32(Y+off);
            cv=vmlaq_f32(cv,av,bv);         
        }
        av=vld1q_f32(tail_c);
        bv=vld1q_f32(tail_d);
        cv=vmlaq_f32(cv,av,bv); 
        r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
        A[n]=vget_lane_f32(vpadd_f32(r,r),0);       
    }
}


#endif
