#pragma once

#include<stdio.h>
#include<stdlib.h>

typedef struct _Image {
  unsigned char*data;
  int W,H;
} Image;

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
