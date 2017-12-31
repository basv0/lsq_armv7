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

#ifndef _MNK_OPTIM_
#define _MNK_OPTIM_

#include <string.h>

inline void mnk_optim(int N, int M, float **X, float **Xt, float *Y, float *A, float **S, float *Q){
    for(int n=0;n<2*M;n++){
        for(int m=n;m<2*M;m++){
            S[n][m]=0.0;
            for(int k=0;k<N;k++){
                S[n][m]+=Xt[n][k]*Xt[m][k];
                if(m!=n) S[m][n]=S[n][m];
            }
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
                for(int k=0;k<4*M;k++){
                    S[j][k]-=ratio*S[i][k];
                }
            }
        }
    }
    for(int i=0;i<2*M;i++){
        a=1.0/S[i][i];
        for(int j=2*M; j<4*M;j++){
            S[i][j]*=a;
        }
    }
    for(int n=0;n<2*M;n++){
        for(int m=0;m<N;m++){
            Q[m]=0.0;
            for(int k=0;k<2*M;k++){
                Q[m] += S[n][k+2*M]*X[m][k];
            }
        }
        A[n]=0;
        for(int m=0;m<N;m++){
            A[n]+=Q[m]*Y[m];
        }
    }
}

#endif
