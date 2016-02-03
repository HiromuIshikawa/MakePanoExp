// gcc TKfilter_3.c image.c matrix.c -lm -O3 -mavx2 -march=native -funroll-loops -fomit-frame-pointer ; ./a.out

#include"image.h"
#include<math.h>
#include<stdio.h>

struct data{
    int x;
    int y;
    struct data *next;
};

#define __rdtsc() ({ long long a,d; asm volatile ("rdtsc":"=a"(a),"=d"(d)); d<<32|a; })
//void ImageFeature(Matrix*im2,Image*im);
struct data *insertion (int x, int y, double add, struct data *head,
                        int n, Matrix* im);
int N=30;

void ImageDrawBox(Image*im,int x,int y){
  int u,v,W=7;
  for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
    IElem(im,x+u,y+v,0)=IElem(im,x+u,y+v,0) + 0xff >> 1;
    IElem(im,x+u,y+v,1)=IElem(im,x+u,y+v,1) >> 1;
    IElem(im,x+u,y+v,2)=IElem(im,x+u,y+v,2) >> 1;
  }
}

void ImageMatrixWrite(Image*im,Matrix*mt,double s){
  int x,y;
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
  //#pragma omp parallel for private (x,v,ix,iy)
  for(y=W+1;y<im->H-W-1;y++) for(x=0;x<im->W-1;x++){
      ixx=iyy=ixy=0;
      for(v=-W;v<=W;v++){
	ix=IElem(im, x+1, y+v,1) - IElem(im, x-1, y+v,1);
	iy=IElem(im, x, y+v+1,1) - IElem(im, x, y+v-1,1);
	ixx+=ix*ix; // ixx だけでなく ixy,iyy も計算する．
	iyy+=iy*iy;
	ixy+=ix*iy;
      }
      DElem(tmp_xx,x,y)=ixx;
      DElem(tmp_yy,x,y)=iyy;
      DElem(tmp_xy,x,y)=ixy;
    }
  //#pragma omp parallel for private (x,u)
  for(y=W+1;y<im->H-W-1;y++) for(x=W+1;x<im->W-W-1;x++){
      ixx=iyy=ixy=0;
      for(u=-W;u<=W;u++){
	//	ix=DElem(tmp, x+u+1, y) - DElem(tmp, x+u-1, y);
	//iy=DElem(tmp, x+u, y+1) - DElem(tmp, x+u, y-1);

	ixx+=DElem(tmp_xx, x+u, y); // ixx だけでなく ixy,iyy も計算する．
	iyy+=DElem(tmp_yy, x+u, y);
	ixy+=DElem(tmp_xy, x+u, y);
      }
      b=ixx+iyy;
      c=ixx*iyy-ixy*ixy;
      DElem(im2,x,y)=(b - sqrt(b*b - 4*c)) / 2; // 実際には [ixx,ixy;ixy,iyy] の小さい方の固有値を入れる．
    }
}

struct data *AssignNewPtr(int x, int y, struct data *next){
    struct data *newptr;
    newptr = (struct data*)malloc(sizeof(struct data));

    if (newptr == NULL) {
        fprintf(stderr,"assign error.\n");
        exit(1);
    }

    newptr->x = x;
    newptr->y = y;
    newptr->next = next;
    return newptr;

}

struct data *MatrixLocalMax(struct data *head, Matrix*im2){


    int x,y,u,v,W=7,n=0;
  for(y=W+1;y<im2->H-W-1;y++) for(x=W+1;x<im2->W-W-1;x++){
      double max=-1;
      for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
	  // (x,y) を中心とする 15x15 の矩形領域内で DElem(im2,x+u,y+v) の最大値を探す．
	  if(max < DElem(im2,x+u,y+v))max = DElem(im2,x+u,y+v);
	}
      if(max == DElem(im2,x,y)) {
	// 最大値が DElem(im2,x,y) と等しいなら，(x,y) を特徴点として記録する．
         head = insertion(x,y,max,head,n,im2);
          if(n < N) n++;
      }
    }
  return head; // 記録した点の数
}

struct data *insertion (int x, int y, double add,struct data *head, int n, Matrix* im){

    struct data *now, *prev;
        if(head == NULL) {/*初めてのデータ登録*/
            return AssignNewPtr(x,y,NULL);
    }

    if (add > DElem(im,head->x,head->y)) {/*先頭に挿入するとき*/
        return AssignNewPtr(x,y,head);
    }

    prev = NULL;
    now = head;

  while (now != NULL) {
      if (add > DElem(im,now->x,now->y)) {/*途中に挿入するとき*/
          prev->next = AssignNewPtr(x,y,prev->next);
          return head;
      }
      prev = now;
      now = now->next;
  }
    if (n <= N) {
        prev->next = AssignNewPtr(x,y,NULL);/*最後に挿入するとき*/
    }
    return head;

}



void PrintFeature (struct data *w, Matrix *im, int fea[][2]){
  int i;
  for (i = 0;i<N && w != NULL;i++, w = w->next) {
    fea[i][0] = w->x;
    fea[i][1] = w->y;
  }
}

void TKfilter(int fea[][2], char *imName){
  Image *im;
  Matrix*im2;

  struct data *kk = NULL;

  im=ImageRead(imName);
  //im3=ImageRead4jpg(av[1]);
  im2=MatrixAlloc(im->H,im->W);


  ImageFeature(im2,im);

  //long long start=__rdtsc();
  kk=MatrixLocalMax(kk,im2);
    //printf("%f msec\n",(__rdtsc()-start)/3.4e+6);
  PrintFeature(kk,im2,fea);
  //ImageMatrixWrite(im,im2,.001);
  //for(i=0; i < N && kk != NULL; i++,kk = kk->next) //ImageDrawBox(im3,kk->x,kk->y);
  //ImageWrite(av[2],im);
  //ImageWrite(av[3],im3);
}
