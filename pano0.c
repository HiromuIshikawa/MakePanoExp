  // gcc -O -lm pano2.c image.c RANSAC.c matrix.c -o PANO && ./a.out && display out.ppm

  #include "image.h"
  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>
  #include <string.h>
  // Following functions must be defined in image.c:
  //   ImageRead, ImageAlloc, ImageClear.
  // In addition, check if your ImageRead can read a jpg file.
  void MatrixSimeqLr(Matrix*mtB,Matrix*mtR);
  Matrix* makeTranceMT(Image *im, Image* im2);
  #define ImageImageProjectionAlpha ImageImageProjectionAlphaNaive


  #define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
  #define Row(_a,_b)     ((_a)->data+(_a)->W*(_b))
  void ImageImageProjectionAlpha(Image*id,Image*is,double a[3][3],double alpha){
    int x,y,u,v;
    double r;
    for(y=0;y<id->H;y++) for(x=0;x<id->W;x++){
      r = 1 / (a[2][0]*x+a[2][1]*y+a[2][2]);
      u = r * (a[0][0]*x+a[0][1]*y+a[0][2]);
      v = r * (a[1][0]*x+a[1][1]*y+a[1][2]);
      if( isInsideImage(is,u,v) ){
        id->data[(y*id->W+x)*3+0] += is->data[(v*is->W+u)*3+0]*alpha,
        id->data[(y*id->W+x)*3+1] += is->data[(v*is->W+u)*3+1]*alpha,
        id->data[(y*id->W+x)*3+2] += is->data[(v*is->W+u)*3+2]*alpha;
      }
    }
  }

  void mult33(double a[3][3], double b[3][3], double c[3][3]){
      int row,col;
      for(row=0;row<3;row++) for(col=0;col<3;col++){
  	a[row][col]=b[row][0]*c[0][col]+
  	            b[row][1]*c[1][col]+
  	            b[row][2]*c[2][col];
      }
  }


  // cc lsq.c -lm -O ; ./a.out
  //#include"mat33.h"

  // basic vector operations


  void MatrixFree(Matrix*mt)
  {
      free(mt->data);
      free(mt);
  }



  void MatrixQRDecompColMajor(Matrix*mtR,Matrix*mt){
    // gram-schmidt orthonormalization (Lt and Q)
    double t, *aT[]= { Row(mt,0),Row(mt,1),Row(mt,2),Row(mt,3),Row(mt,4),Row(mt,5),Row(mt,6),Row(mt,7)};//*aT[mt->W];
    int W = mt->W;
    int i,j;
    MatrixClear(mtR);
    for (i = 0; i < 8;i++) {
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

  Matrix* lsq(int n, double xy [][2], double uv [][2]){
    Matrix *cmA, *vt, *mtR, *tmp;
    int i;
    double z=1;
  //0&1 first
  /* double xy[][2]={ // from 0.jpg
      256,218,
      347,220,
      263,367,
      413,315,
      337,518 //5th point
    },uv[][2]={ // from 1.jpg
      371,230,
      463,230,
      383,379,
      530,327,
      454,530 //5th point
      };*/
    cmA=MatrixAlloc(8,2*n);
    vt=MatrixAlloc(1,2*n);

    // create A (col-major)
    for(i=0;i<n;i++){
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
    MatrixFree(cmA);
    MatrixFree(vt);
    MatrixFree(mtR);
      return tmp;
  }


  int main(int ac, char *av[]){
    Image *im,*im2,*im11,*im22;
    Matrix *ans;
    int i,j;
    im2=ImageAlloc(1024,768);
    ImageClear(im2);
    {
      double m0d[][3]={
        1,0,-100,
        0,1,-100,
        0,0,1
      };
      im=ImageRead4jpg(av[1]);
      ImageImageProjectionAlpha(im2,im,m0d,.5);
    }

    {
        double m10[3][3], m1d[3][3];
      double m0d[][3]={
        1,0,-100,
        0,1,-100,
        0,0,1
      };
      im11=ImageRead4jpg(av[1]);
      im22=ImageRead4jpg(av[2]);
      ans=MatrixAlloc(1,8);
      ans=makeTranceMT(im11, im22);
        for(i=0;i<ans->H;i++){
  	for(j=0;j<ans->W;j++)
  	  m10[i][j] = Elem(ans,i,j);
        }
  	m10[2][2] = 1;
      mult33(m1d,m10,m0d);
      im=ImageRead4jpg(av[2]);
      ImageImageProjectionAlpha(im2,im,m1d,.5);
    }

    ImageWrite(av[3],im2);

    return 0;
  }
