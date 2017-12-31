#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>

const int Vs=4;

inline void mnk(int N, int M, int M2, 
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

int main(int argc, char *argv[]){
    if(argc<3){	
        printf("Использование $%s M N\n",argv[0]);
        printf("M-порядок, N-колич.точек\n");
        return 0;
    }
    M=atoi(argv[1]);
    N=atoi(argv[2]);
    int M2 = (2*M)/Vs+(int)((2*M)%Vs!=0);
    
    X=new float*[N];
    Xt=new float*[2*M];
    Xv=new float32x4_t*[N];
    A=new float[M2*Vs];
    Y=new float[N];
    
    Sv=new float32x4_t*[M2];
    Qv=new float32x4_t[N];
    
    for(int i=0;i<M2;i++){
        Sv[i]=new float32x4_t[4*M];
    }
    
    for(int i=0;i<N;i++){
        X[i]= new float[2*M];
        Xv[i]=new float32x4_t[M2];
    }
    
    for(int j=0;j<2*M;j++){
        Xt[j]=new float[N];
    }

    for(int i=0;i<N;i++) 
        Y[i]=  (rand()%5000000L-2500000L)/1000.0;//yV[i+M];

    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            X[i][j]=  (rand()%5000000L-2500000L)/1000.0;//(-1.0)*yV[i+M-j-1];
            Xt[j][i]=X[i][j];
        }
    }

    for(int i=0; i<N;i++){
        for(int j=0;j<M;j++){
            X[i][j+M]=  (rand()%5000000L-2500000L)/1000.0;//uV[i+M-j-1];
            Xt[j+M][i]=X[i][j+M];
        }
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

    mnk(N, M, M2, X, Xt, Y, A, Xv, Sv, Qv);
    
    clock_gettime(CLOCK_REALTIME, &t2);
    dt = (t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_nsec-t1.tv_nsec)/1000.0;
    for(int i=0;i<2*M;i++){
        fprintf(stderr,"%f\n",A[i]);
    }
    fprintf(stderr,"t=%f\n",dt);
        printf("%d %d %f\n",M,N,dt);
    
    for(int i=0;i<N;i++){
        delete[] X[i];
        delete[] Xv[i];
    }
    
    for(int j=0;j<2*M;j++){
        delete[] Xt[j];
    }
    for(int i=0;i<M2;i++){
        delete[] Sv[i];
    }
    
    delete[] X;
    delete[] Xv;
    delete[] Xt;
    delete[] A;
    delete[] Y;
    delete[] Qv;
    delete[] Sv;
    
    return 0;
}
