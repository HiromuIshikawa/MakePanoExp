#pragma once

#include<stdio.h>
#include<stdlib.h>

typedef struct _Image {
  unsigned char*data;
  int W,H;
} Image;

typedef struct {
  double *data;
  int W,H;
} Matrix;

Image*ImageAlloc(int _W,int _H);
Image*ImageRead(const char*name);
Image*ImageRead4jpg(const char*name);
void ImageFree(Image*im);
void ImageWrite(const char*name,Image*im);
void ImageWrite2jpg(const char*name, Image*im);
void ImageClear(Image*im);
#define IElem(_im,_x,_y,_c) (_im) -> data[(_y)*(_im)->W*3 +(_x)*3+(_c)]
#define isInsideImage(is,u,v) (((unsigned)u<is->W)&&((unsigned)v<is->H))
void ImageImageProjectionAlpha(Image*id,Image*is,double a[3][3],double alpha);
void ImageDrawBox(Image*im,int x,int y);


#define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
#define Row(_a,_b)     ((_a)->data+(_a)->W*(_b))
#define DElem(_a,_b,_c)  (_a)->data[(_a)->W*(_c)+(_b)]

     double VP(double*a,double*b,int N);
     void VSS(double*d,double s,int N);
     void VSA(double*d,double*a,double s,int N);
     Matrix*MatrixAlloc(int _H,int _W);
     void MatrixClear(Matrix*mt);
     void MatrixCopy(Matrix*mtD,Matrix*mt);
     void MatrixCopyT(Matrix*mtD,Matrix*mt);
     void MatrixPrint(Matrix*mt);
     void MatrixMultT(Matrix*mtD,Matrix*mtA,Matrix*mtB);
