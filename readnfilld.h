/*
macro takes data from text file having experimental points
(time/input/output) describing dynamic object, allocates memory
for the matrix of experiment X, its transpose counterpart Xt,
vector Y and fills them with the data.
*/
#ifndef _READNFILL_
#define _READNFILL_

    #define readNfill                                \
    M=atoi(argv[1]);                                 \
    yu = fopen(argv[2],"r");                         \
    if(!yu){                                         \
        printf("file %s not found\n", argv[2]);      \
        return 2;                                    \
    }                                                \
    std::vector<double> yV,uV,tV;                    \
    char syu[80];                                    \
    double u,y,t;                                    \
    while(fgets(syu,80,yu)){                         \
        sscanf(syu,"%lf %lf %lf",&t,&u,&y);          \
        yV.push_back(y);                             \
        uV.push_back(u);                             \
        tV.push_back(t);                             \
    }                                                \
    N=(yV.size()-M);                                 \
    X=new double*[N];                                \
    Xt=new double*[2*M];                             \
    A=new double[2*M];                               \
    Y=new double[N];                                 \
    for(int i=0;i<N;i++){                            \
        X[i]= new double[2*M];                       \
    }                                                \
    for(int j=0;j<2*M;j++){                          \
        Xt[j]=new double[N];                         \
    }                                                \
    for(int i=0;i<N;i++){                            \
        Y[i]= yV[i+M];                               \
    }                                                \
    for(int i=0;i<N;i++){                            \
        for(int j=0;j<M;j++){                        \
            X[i][j]=(-1.0)*yV[i+M-j-1];              \
            Xt[j][i]=X[i][j];                        \
        }                                            \
    }                                                \
    for(int i=0; i<N;i++){                           \
        for(int j=0;j<M;j++){                        \
            X[i][j+M]=uV[i+M-j-1];                   \
            Xt[j+M][i]=X[i][j+M];                    \
        }                                            \
    }                                                \
    yV.clear();                                      \
    uV.clear();                                      \
    tV.clear();                                      \
    
#endif
