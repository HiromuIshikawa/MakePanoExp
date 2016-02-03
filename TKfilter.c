// gcc TKfilter.c image.c matrix.c -lm -O3 -mavx2 -march=native -funroll-loops -fomit-frame-pointer ; ./a.out

#include"image.h"
#include<math.h>

#define __rdtsc() ({ long long a,d; asm volatile ("rdtsc":"=a"(a),"=d"(d)); d<<32|a; })
#define DElem(_a,_b,_c)  (_a)->data[(_a)->W*(_c)+(_b)]
//void ImageFeature(Matrix*im2,Image*im);

void ImageDrawBox(Image*im,int x,int y){
  int u,v,W=7;
  for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
    IElem(im,x+u,y+v,0)=IElem(im,x+u,y+v,0) + 0xff >> 1;
    IElem(im,x+u,y+v,1)=IElem(im,x+u,y+v,1) >> 1;
    IElem(im,x+u,y+v,2)=IElem(im,x+u,y+v,2) >> 1;
  }
}

void ImageMatrixWrite(Image*im,Matrix*mt,double s){
  int x,y,p;
  for(y=0;y<mt->H;y++) for(x=0;x<mt->W;x++){
    double tmp=mt->data[y*mt->W+x]*s;
    if(tmp>255) tmp=255;
    im->data[(y*im->W+x)*3+0]=
    im->data[(y*im->W+x)*3+1]=
    im->data[(y*im->W+x)*3+2]=tmp;
  }
}

void ImageFeature(Matrix*im2,Image*im){
  int x,y,u,v,W=7,ix,iy;
  double ixx,ixy,iyy,b,c;
  Matrix*tmp_xx, *tmp_yy, *tmp_xy;
  tmp_xx=MatrixAlloc(im->H,im->W);
  tmp_yy=MatrixAlloc(im->H,im->W);
  tmp_xy=MatrixAlloc(im->H,im->W);
    for(y=W+1;y<im->H-W-1;y++) for(x=W+1;x<im->W-W-1;x++){
      ixx=iyy=ixy=0;
      for(v=-W;v<=W;v++){
	ix=IElem(im, x+1, y+v,1) - IElem(im, x-1, y+v,1);
	iy=IElem(im, x, y+v+1,1) - IElem(im, x, y+v-1,1);
	ixx+=ix*ix;
	iyy+=iy*iy;
	ixy+=ix*iy;
      }
      DElem(tmp_xx,x,y)=ixx;
      DElem(tmp_yy,x,y)=iyy;
      DElem(tmp_xy,x,y)=ixy;
    }
  for(y=W+1;y<im->H-W-1;y++) for(x=W+1;x<im->W-W-1;x++){
      ixx=iyy=ixy=0;
      for(u=-W;u<=W;u++){
	
	ixx+=DElem(tmp_xx, x+u, y);
	iyy+=DElem(tmp_yy, x+u, y);
	ixy+=DElem(tmp_xy, x+u, y);
        }
      b=ixx+iyy;
      c=ixx*iyy-ixy*ixy;
      DElem(im2,x,y)=(b - sqrt(b*b - 4*c)) / 2;
  }
}

int MatrixLocalMax(int w[][2], Matrix*im2){
  int x,y,u,v,W=7,n=0,a;
  for(y=W+1;y<im2->H-W-1;y++) for(x=W+1;x<im2->W-W-1;x++){
    double max=-1;
    for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
      // (x,y) を中心とする 15x15 の矩形領域内で DElem(im2,x+u,y+v) の最大値を探す．
    }
    // 最大値が DElem(im2,x,y) と等しいなら，(x,y) を特徴点として記録する．
    // a=n++; w[a][0]=x; w[a][1]=y;
  }
  return n; // 記録した点の数
}


int main(int ac,char**av){
  Image *im;
  Matrix*im2;
  int kk[9999][2], kw,i;

  if(ac<2) return 1;

  im=ImageRead4jpg(av[1]);
  im2=MatrixAlloc(im->H,im->W);
long long start=__rdtsc();
  ImageFeature(im2,im);
printf("%f msec\n",(__rdtsc()-start)/3.4e+6);
  kw=MatrixLocalMax(kk,im2);

  ImageMatrixWrite(im,im2,.001);
  for(i=0;i<kw;i++) ImageDrawBox(im,kk[i][0],kk[i][1]);
  ImageWrite("out3.ppm",im);
}
