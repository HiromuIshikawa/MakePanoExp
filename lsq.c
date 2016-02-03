// cc lsq.c -lm -O ; ./a.out

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
//#include"mat33.h"

// basic vector operations
double VP(double*a,double*b,int N){
  double s=0;
  int i;
  for(i=0;i<N;i++) s += a[i] * b[i] ;
  return s;
}

void VSS(double*d,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] *= s;
}

void VSA(double*d,double*a,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] += a[i] * s;
}

// matrix structure and functions
typedef struct {
  double *data;
  int W,H;
} Matrix;

#define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
#define Row(_a,_b)     ((_a)->data+(_a)->W*(_b))

Matrix*MatrixAlloc(int _H,int _W){
  Matrix*mt=(Matrix*)malloc(sizeof(Matrix));;
  mt->W = _W,
  mt->H = _H;
  mt->data=(double*)malloc(mt->W*mt->H*sizeof(double));
  return mt;
}

void MatrixClear(Matrix*mt){
  memset(mt->data,0,mt->W*mt->H*sizeof(double));
}


void MatrixCopy(Matrix*mtD,Matrix*mt){
  memmove(mtD->data,mt->data,mt->W*mt->H*sizeof(double));
}

void MatrixCopyT(Matrix*mtD,Matrix*mt){
  int i,j;
  for(i=0;i<mtD->H;i++)
    for(j=0;j<mtD->W;j++)
      Elem(mtD,i,j) = Elem(mt,j,i);
}

void MatrixPrint(Matrix*mt){
  int i,j;
  for(i=0;i<mt->H;i++){
    for(j=0;j<mt->W;j++)
      printf("%f ",Elem(mt,i,j));
    printf("\n");
  }
}

void MatrixMultT(Matrix*mtD,Matrix*mtA,Matrix*mtB){
  // D = A B^T
  int i,j;
  for(i=0;i<mtA->H;i++)
    for(j=0;j<mtB->H;j++)
      Elem(mtD,i,j) = VP( Row(mtA,i), Row(mtB,j), mtA->W);
}


void MatrixQRDecompColMajor(Matrix*mtR,Matrix*mt){
  // gram-schmidt orthonormalization (Lt and Q)
  double t, *aT[mt->W];
  int W = mt->W;
  int i,j;
  MatrixClear(mtR); 
  for (i = 0; i < W;i++) {
    aT[i] = Row(mt,i);
    for (j = 0; j < i; j++) {
      Elem(mtR,j,i) = t = VP(aT[j], aT[i], W);
      VSA(aT[i], aT[j], -t, W);
    }
    Elem(mtR,i,i) = t = sqrt(VP(aT[i],aT[i],W));
    VSS(aT[i], 1/t, W); 
  }




  /*

  Elem(mtR,0,0) = t = sqrt(VP(aT[0],aT[0],W));
  VSS(aT[0], 1/t, W);


  Elem(mtR,0,1) = t = VP(aT[0], aT[1], W);
  VSA(aT[1], aT[0], -t, W);

  Elem(mtR,1,1) = t = sqrt(VP(aT[1],aT[1],W));
  VSS(aT[1], 1/t, W);

///////////
  Elem(mtR,0,2) = t = VP(aT[0], aT[2], W);
  VSA(aT[2], aT[0], -t, W);

  Elem(mtR,1,2) = t = VP(aT[1], aT[2], W);
  VSA(aT[2], aT[1], -t, W);

  Elem(mtR,2,2) = t = sqrt(VP(aT[2],aT[2],W));
  VSS(aT[2], 1/t, W);

////// 以下略
*/
}

void MatrixSimeqLr(Matrix*mtB,Matrix*mtR){
  // B = B L^{-1}
  double * B = Row(mtB,0);
  int i,j;
  for (i = mtB->W - 1; i >= 0; i--) {
    for (j = i+1; j < mtB->W; j++) {
      B[i] = B[i]-B[j]*Elem(mtR,i,j);
    }
    B[i] = B[i] / Elem(mtR,i,i);
  }

  /*
  B[7] =  B[7] / Elem(mtR,7,7);
  B[6] = (B[6]-B[7]*Elem(mtR,6,7)) / Elem(mtR,6,6);
  B[5] = (B[5]-B[6]*Elem(mtR,5,6)-B[7]*Elem(mtR,5,7)) / Elem(mtR,5,5);
///// 以下略
*/
  
}


main(){
  Matrix *cmA, *vt, *mtR, *tmp;
  int i;
  double z=1;

  double xy[][2]={ // from 0.jpg
    148,537,
    347,220,
    263,367,
    413,315,
  },uv[][2]={ // from 1.jpg
    371,230,
    463,230,
    383,379,
    530,327,
  };
  double ans[3][3];
  cmA=MatrixAlloc(8,8);
  vt=MatrixAlloc(1,8);

  // create A (col-major)
  for(i=0;i<4;i++){
    cmA->data[cmA->W*0+(i*2  )]=z*xy[i][0];
    cmA->data[cmA->W*1+(i*2  )]=z*xy[i][1];
    cmA->data[cmA->W*2+(i*2  )]=z*z;
    cmA->data[cmA->W*3+(i*2  )]=0;
    cmA->data[cmA->W*4+(i*2  )]=0;
    cmA->data[cmA->W*5+(i*2  )]=0;
    cmA->data[cmA->W*6+(i*2  )]=-xy[i][0]*uv[i][0];
    cmA->data[cmA->W*7+(i*2  )]=-xy[i][1]*uv[i][0];
    cmA->data[cmA->W*0+(i*2+1)]=0;
    cmA->data[cmA->W*1+(i*2+1)]=0;
    cmA->data[cmA->W*2+(i*2+1)]=0;
    cmA->data[cmA->W*3+(i*2+1)]=z*xy[i][0];
    cmA->data[cmA->W*4+(i*2+1)]=z*xy[i][1];
    cmA->data[cmA->W*5+(i*2+1)]=z*z;
    cmA->data[cmA->W*6+(i*2+1)]=-xy[i][0]*uv[i][1];
    cmA->data[cmA->W*7+(i*2+1)]=-xy[i][1]*uv[i][1];
    vt->data[i*2  ]=z*uv[i][0];
    vt->data[i*2+1]=z*uv[i][1];
  }

  // solve Least-squares equation
  mtR=MatrixAlloc(8,8);
  MatrixQRDecompColMajor(mtR,cmA);
  tmp=MatrixAlloc(1,8);
  MatrixMultT(tmp,vt,cmA);
  MatrixSimeqLr(tmp,mtR);
  MatrixPrint(tmp);
}
