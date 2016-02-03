#include <string.h>
#include"image.h"

Image*ImageAlloc(int W,int H){
  Image*im=(Image*)malloc(sizeof(Image));
  im->W=W;
  im->H=H;
  im->data=(unsigned char*)malloc(W*H*3);
  return im;
}

Image*ImageRead(const char*name){
  int W,H;
  char *imtype;
  Image *im;
  FILE *fp;
  imtype = strstr(name, ".");

  if (strcmp(imtype, ".ppm") == 0) {
    if ((fp=fopen(name,"rb")) == NULL) {
      fprintf(stderr, "%s cannot open.\n", name);
    }
    fscanf(fp,"%*s%d%d%*s%*c",&W,&H);
    im=ImageAlloc(W,H);
    fread(im->data,1,W*H*3,fp);
    fclose(fp);
  } else if (strcmp(imtype, ".jpg") == 0) {
    //printf("%s type is %s\n", name, imtype);
    char cmd[strlen(name)+6];
    sprintf(cmd, "djpeg %s",name);
    if ((fp=popen(cmd,"r")) == NULL) {
      fprintf(stderr, "%s cannot open.\n", name);
    }
    fscanf(fp,"%*s%d%d%*s%*c",&W,&H);
    im=ImageAlloc(W,H);
    fread(im->data,1,W*H*3,fp);
    pclose(fp);
  } else {
    fprintf(stderr, "this imagetype is not availablable\n");
    im = NULL;
  }

  //im=ImageAlloc(W,H);
  return im;
}
Image*ImageRead4jpg(const char*name){//画像ファイルを読み込む関数
  int W,H;
  Image*im;
  char cmd[strlen(name)+6];
  sprintf(cmd, "djpeg %s",name);
  FILE*fp=popen(cmd,"r");
  fscanf(fp,"%*s%d%d%*s%*c",&W,&H);//`*'を指定すると変数に書き込まれない
  im=ImageAlloc(W,H);
  fread(im->data,1,W*H*3,fp);
  pclose(fp);
  return im;
}

void ImageWrite(const char*name,Image*im)
{
  FILE *fp;
  char *imtype;
  imtype = strstr(name, ".");
  if (strcmp(imtype, ".ppm") == 0) {
    if ((fp = fopen(name,"wb")) == NULL) {
      fprintf(stderr, "%s cannot open.\n", name);
    }else{
      fprintf(fp, "P6 %d %d 255\n",im->W, im->H);
      fwrite(im->data,1, im->W*im->H*3, fp);
      fclose(fp);
    }
  } else if (strcmp(imtype, ".jpg") == 0) {
    //fprintf(stderr,"%s type is %s\n", name, imtype);
    char cmd[strlen(name)+8];
    sprintf(cmd, "cjpeg > %s",name);
    if ((fp = popen(cmd,"w")) == NULL) {
      fprintf(stderr, "%s cannot open.\n", name);
    } else {
      fprintf(fp, "P6 %d %d 255\n",im->W, im->H);
      fwrite(im->data,1, im->W*im->H*3, fp);
      pclose(fp);
    }
  } else {
    fprintf(stderr, "this imagetype is not availablable\n");
  }


}

void ImageWrite2jpg(const char*name,Image*im)/* jpg 出力関数*/
{
  FILE *fp;
  char cmd[strlen(name)+8];
  sprintf(cmd, "cjpeg > %s",name);
  fp = popen(cmd,"w");
  fprintf(fp, "P6 %d %d 255\n",im->W, im->H);
  fwrite(im->data,1, im->W*im->H*3, fp);
  pclose(fp);

}

void ImageClear(Image*im){
  int x,y;
  for(y=0;y<im->H;y++) for(x=0;x<im->W;x++){
      im->data[(y*im->W+x)*3+0] = 0;
      im->data[(y*im->W+x)*3+1] = 0;
      im->data[(y*im->W+x)*3+2] = 0;
    }

}
