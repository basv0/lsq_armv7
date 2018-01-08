#ifndef _FREEMEM_
#define _FREEMEM_

#define freeMem                 \
    for(int i=0;i<N;i++){       \
	delete[] X[i];          \
    }                           \
    for(int j=0;j<2*M;j++){     \
	delete[] Xt[j];         \
    }                           \
    delete[] X;                 \
    delete[] Xt;                \
    delete[] A;                 \
    delete[] Y;                 \

#endif
