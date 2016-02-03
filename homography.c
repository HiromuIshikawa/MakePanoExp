// gcc -O homography.c image.c && ./a.out 1.jpg out.ppm && display out.ppm

#include"image.h"

// Following functions must be defined in image.c:
//   ImageRead, ImageAlloc, ImageClear, ImageWrite.
// In addition, check whether your ImageRead can read a jpg file.

void ImageImageProjection(Image*id,Image*is,double a[3][3]){
  int x,y,u,v;
  double r;
  for(y=0;y<id->H;y++) for(x=0;x<id->W;x++){
    r = 1 / (a[2][0]*x+a[2][1]*y+a[2][2]);
    u = r * (a[0][0]*x+a[0][1]*y+a[0][2]);
    v = r * (a[1][0]*x+a[1][1]*y+a[1][2]);
    if( isInsideImage(is,u,v) ){
      id->data[(y*id->W+x)*3+0] = is->data[(v*is->W+u)*3+0],
      id->data[(y*id->W+x)*3+1] = is->data[(v*is->W+u)*3+1],
      id->data[(y*id->W+x)*3+2] = is->data[(v*is->W+u)*3+2];
    }
  }
}

int main(int argc,char*argv[]){
  Image *im,*im2;
  double a[][3]={
     .866 , -.5  , 160,
     .5   , .866 , -300,
     0    , 0    , 1
  };

  double b[][3]={
     .866 , -.5  ,  160,
     .5   , .866 , -300,
     -0.001, 0   , 1
  };
  double c[][3]={
     .0 , -1.0  ,  770,
     1.0   , .0 , -300,
     0, 0   , 1
  };
  im=ImageRead4jpg(argv[1]);

  im2=ImageAlloc(1024,768);
  ImageClear(im2);
  ImageImageProjection(im2,im,c);
  ImageWrite(argv[2],im2);

  return 0;
}
