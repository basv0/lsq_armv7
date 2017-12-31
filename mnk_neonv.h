/**************************************************************************
*   Copyright (C) 2017 by Basil Olonichev                                 *
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

#ifndef _NEONV_
#define _NEONV_

#include <string.h>
#include <arm_neon.h>

const int Vs=128/sizeof(float)/8;

inline void mnk_neonv(int N, int M, int M2, 
                float **X, float**Xt, float *Y, float* A,
                float32x4_t **Xv, float32x4_t **Sv, float32x4_t *Qv){
    float zero=0.0;
    float32x4_t av;
    float32x4_t bv;    
    float fm[Vs];
    
    for(int m=0;m<M2;m++){
        for(int j=0;j<2*M;j++){
            Sv[m][j]=vld1q_dup_f32(&zero);
            for(int n=0;n<N;n++){
                av=Xv[n][m];
                bv=vmulq_n_f32(av,Xt[j][n]);
                Sv[m][j]=vaddq_f32(Sv[m][j],bv);
            }
        }
    }
    
    for(int i=0;i<2*M;i++){
        int m1=i/Vs;
        int m2=i%Vs;
        vst1q_f32(fm,Sv[m1][i+2*M]);
        fm[m2]=1.0;
        Sv[m1][i+2*M]=vld1q_f32(fm);
    }
    
    for(int i=0;i<2*M;i++){
        int m1=i/Vs;
        int m2=i%Vs;
        vst1q_f32(fm,Sv[m1][i]);
        float z=1.0/fm[m2];
        for(int j=0;j<2*M;j++){
            if(i!=j){
                int m3=j/Vs;
                int m4=j%Vs;
                vst1q_f32(fm,Sv[m3][i]);
                float ratio=fm[m4]*z;
                for(int k=0;k<4*M;k++){
                    float fi[4];
                    float fj[4];
                    vst1q_f32(fj,Sv[m3][k]);
                    vst1q_f32(fi,Sv[m1][k]);
                    fj[m4]=fj[m4]-fi[m2]*ratio;
                    Sv[m3][k]=vld1q_f32(fj);
                }
            }
        }
    }
    for(int i=0;i<2*M;i++){
        int m1=i/Vs;
        int m2=i%Vs;
        vst1q_f32(fm,Sv[m1][i]);
        float a=1.0/fm[m2];
        for(int j=2*M;j<4*M;j++){
            vst1q_f32(fm,Sv[m1][j]);
            fm[m2]*=a;
            Sv[m1][j]=vld1q_f32(fm);
        }
    }
    for(int m=0;m<M2;m++){
        av=vld1q_dup_f32(&zero);
        for(int n=0;n<N;n++){
            Qv[n]=vld1q_dup_f32(&zero);
            for(int i=0;i<2*M;i++){
                bv=vmulq_n_f32(Sv[m][i+2*M],X[n][i]);
                Qv[n]=vaddq_f32(Qv[n],bv);
            }
        }
        for(int n=0;n<N;n++){
            bv = vmulq_n_f32(Qv[n],Y[n]);
            av = vaddq_f32(av,bv);
        }
        vst1q_f32(A+m*Vs,av);
    }
}

#endif
