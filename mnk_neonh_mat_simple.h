/**************************************************************************
*   Copyright (C) 2019 by Basil Olonichev                                 *
*   basv0255@gmail.com                                                    *
*   v_olonichev@ksu.edu.ru
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

const int Vs=128/sizeof(float)/8; //size of vector

//left-overs are processed as scalar values
inline void mnk(int N, int M, float **X, float **Xt, float *Y, float *A, float **S, float *Q){
    int D=N/Vs;
    float zero=0.0;
    float32x4_t av;
    float32x4_t bv;
    float32x4_t cv;
    float32x2_t r;
    for(int n=0;n<2*M;n++){
        for(int m=n;m<2*M;m++){
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<D;k++){
                int off=k*Vs;
                av=vld1q_f32(Xt[n]+off);
                bv=vld1q_f32(Xt[m]+off);
                cv=vmlaq_f32(cv,av,bv);         
            }
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            S[n][m]=vget_lane_f32(vpadd_f32(r,r),0);
            for(int k=D*Vs;k<N;k++){
                S[n][m]+=Xt[n][k]*Xt[m][k];
            }
            S[m][n]=S[n][m];
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
    for(int n=0;n<2*M;n++){
        for(int m=0;m<N;m++){
            cv=vld1q_dup_f32(&zero);
            for(int k=0;k<D;k++){
                int off=k*Vs;
                av=vld1q_f32(S[n]+2*M+off);
                bv=vld1q_f32(X[m]+off);
                cv=vmlaq_f32(cv,av,bv);         
            }
            r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
            Q[m]=vget_lane_f32(vpadd_f32(r,r),0);       
            for(int k=D*Vs;k<2*M;k++){
                Q[m]+=S[n][2*M+k]*X[m][k];
            }
        }
        int D1=N/Vs;
        cv=vld1q_dup_f32(&zero);
        for(int m=0;m<D1;m++){
            int off=m*Vs;
            av=vld1q_f32(Q+off);
            bv=vld1q_f32(Y+off);
            cv=vmlaq_f32(cv,av,bv);         
        }
        r=vadd_f32(vget_high_f32(cv),vget_low_f32(cv));
        A[n]=vget_lane_f32(vpadd_f32(r,r),0);
        for(int m=D1*Vs;m<N;m++){
            A[n]+=Q[m]*Y[m];
        }
    }
}
#endif
