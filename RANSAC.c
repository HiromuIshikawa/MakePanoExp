  /*

  gcc -O RANSAC.c image.c matrix.c -lm -o greedy; ./a.out 0.jpg 2.jpg
  */
  #include"image.h"
  #define INFINITY (1/0.)
  Matrix* lsq(int n, double xy [][2], double uv [][2]);
  void MatrixFree(Matrix*mt);


  int ransacMethod(int N, int ansAry[][4],int m[][2], int x1[][2],int x2[][2]) {
    int i,j,roop,snm=0,fnm,rndAry[N];
    Matrix *ans;
    double tmp[3][3],w[100][4];

    for (i = 0; i < N; i++) rndAry[i] = i;

    for (roop = 0; roop <= 1000; roop++) {
      fnm = 0;
      for(i=0;i < 4;i++){
        int t;
        j=(int)(drand48()*N);
        t=rndAry[i]; rndAry[i]=rndAry[j]; rndAry[j]=t;
      }
      double xy[][2]={ // from the first image
        x1[m[rndAry[0]][0]][0], x1[m[rndAry[0]][0]][1],
        x1[m[rndAry[1]][0]][0], x1[m[rndAry[1]][0]][1],
        x1[m[rndAry[2]][0]][0], x1[m[rndAry[2]][0]][1],
        x1[m[rndAry[3]][0]][0], x1[m[rndAry[3]][0]][1]
      },uv[][2]={ // from the second image
        x2[m[rndAry[0]][1]][0], x2[m[rndAry[0]][1]][1],
        x2[m[rndAry[1]][1]][0], x2[m[rndAry[1]][1]][1],
        x2[m[rndAry[2]][1]][0], x2[m[rndAry[2]][1]][1],
        x2[m[rndAry[3]][1]][0], x2[m[rndAry[3]][1]][1]
      };
      ans = MatrixAlloc(1,8);
      ans = lsq(4, xy, uv);
      for(i=0;i<ans->H;i++){
        for(j=0;j<ans->W;j++){
  	tmp[i][j] = Elem(ans,i,j);
        }
      }
      tmp[2][2] = 1;
      MatrixFree(ans);
      for(i=0;i<30;i++){
        double x=x1[m[i][0]][0], y=x1[m[i][0]][1],
  	u=x2[m[i][1]][0], v=x2[m[i][1]][1],
  	dx,dy,r,xu,yv;
        r = 1 / (tmp[2][0]*x+tmp[2][1]*y+tmp[2][2]);
        xu = r * (tmp[0][0]*x+tmp[0][1]*y+tmp[0][2]);
        yv = r * (tmp[1][0]*x+tmp[1][1]*y+tmp[1][2]);
        // dx, dy は (x,y) を変換した座標と (u,v) の差．
        dx = xu - u;
        dy = yv - v;
        if(dx*dx+dy*dy < 5){
  	w[fnm][0] = x;
  	w[fnm][1] = y;
  	w[fnm][2] = u;
  	w[fnm][3] = v;
  	fnm++;
        }
      }
      if (fnm > snm) {
        snm = fnm;
        for(i = 0; i < snm; i++) {
  	ansAry[i][0] = w[i][0];
  	ansAry[i][1] = w[i][1];
  	ansAry[i][2] = w[i][2];
  	ansAry[i][3] = w[i][3];
        }
      }
    }
    return snm;
  }


  int greedyMethod(int match[][2],Matrix*mt,
  		 Image*im ,int x1[][2],int N1,
  		 Image*im2,int x2[][2],int N2){
    int i,j,k,ji=0,n=0;

    for(i=0;i<N1;i++){
      double sm=INFINITY,t;
      for(j=0;j<N2;j++){
        t=Elem(mt,i,j);
        // printf("%f ",t);
        if( sm > t ) sm = t, ji=j;
      }
      match[n][0]=i; match[n][1]=ji; n++;
      /*printf("%d,%d,%d,%d,\n",
  	   x1[i][0],x1[i][1],
  	   x2[ji][0],x2[ji][1]);*/
      for(k=0;k<N1;k++) Elem(mt,k,ji) = INFINITY;
    }

    return n;
    }


  int matchMethod2(int match[][2],Matrix*mt,
  		 Image*Im ,int x1[][2],int N1,
  		 Image*Im2,int x2[][2],int N2) {
    int i,j,n,in=0,jn=0,k,cnt;
    if(N1 > N2){
      cnt = N2;
    } else {
      cnt = N1;
    }
    for (n = 0; n < cnt; n++) {
      double sm = INFINITY,t;
      for(i = 0; i < N1; i++) {
        for(j = 0; j < N2; j++) {
  	t = Elem(mt,i,j);
  	if(sm > t){
  	  sm = t;
  	  in = i;jn = j;
  	}
        }
      }
      match[n][0] = in; match[n][1] = jn;
      for(k=0;k<N1;k++) Elem(mt,k,jn) = INFINITY;
      for(k=0;k<N2;k++) Elem(mt,in,k) = INFINITY;
    }
    return n;
  }


  double ImageSSD(Image*im,int x1,int y1, Image*im2,int x2,int y2){
    int i,j,W=7;
    double sr=0,sg=0,sb=0,dr,dg,db;
    for(i=-W;i<=W;i++) for(j=-W;j<=W;j++){
      dr  = IElem(im, x1+j, y1+i, 0) - IElem(im2, x2+j , y2+i, 0);
      dg  = IElem(im, x1+j, y1+i, 1) - IElem(im2, x2+j , y2+i, 1);
      db  = IElem(im, x1+j, y1+i, 2) - IElem(im2, x2+j , y2+i, 2);
      sr += dr*dr;
      sg += dg*dg;
      sb += db*db;
    }
    return sr+sg+sb;
  }


  void calcSSDtable(Matrix*mt,
  		  Image*im ,int x1[][2],int N1,
  		  Image*im2,int x2[][2],int N2){
    int i,j;
    for(i=0;i<N1;i++){
      for(j=0;j<N2;j++){
        Elem(mt,i,j) = ImageSSD(im ,x1[i][0],x1[i][1],
  			      im2,x2[j][0],x2[j][1]);
      }
    }
  }


  Matrix* makeTranceMT(Image *im, Image* im2){
    Matrix *mt, *ans;
    int i,j,nm,snm;
    int match[999][2],ansAry[100][4];
    int x1[][2]={
  #include"0.fea"
    }, N1=30;
    int x2[][2]={
  #include"1.fea"
    }, N2=30;

    mt=MatrixAlloc(N1,N2);
    calcSSDtable(mt,im,x1,N1,im2,x2,N2);

    nm=matchMethod2(match,mt,im,x1,N1,im2,x2,N2); // 特徴点の対応付け
    /* for(i = 0; i < nm; i++)
      printf("%d:(%d,%d,%d,%d)\n",i,x1[match[i][0]][0],x1[match[i][0]][1],
      x2[match[i][1]][0],x2[match[i][1]][1]);*/

    snm = ransacMethod(N1, ansAry, match, x1, x2);
    double xy[snm][2], uv[snm][2];
    for(j = 0; j < snm; j++) {
      xy[j][0] = ansAry[j][0];
      xy[j][1] = ansAry[j][1];
      uv[j][0] = ansAry[j][2];
      uv[j][1] = ansAry[j][3];
      printf("(%f %f %f %f)\n", xy[j][0],xy[j][1],uv[j][0],uv[j][1]);
    }
    ans = MatrixAlloc(1,8);
    printf("-----------------%d\n",snm);
    ans = lsq(snm,xy,uv);
    return ans;
  }
