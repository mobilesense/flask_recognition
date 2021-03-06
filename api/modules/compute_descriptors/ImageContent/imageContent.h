// Definition de la classe ImageContent
#ifndef _imageContent_h_
#define _imageContent_h_

#include <cstdlib>
#include <fstream>
#include <string.h>
#include <cstdio>
#include <cmath>
#include <iostream>
//#include "/usr/local/libpng/include/png.h"
#include <png.h>

using namespace std;

#define COL 5
#define GRAY 1
#define FLOAT 2
#define UCHAR 1
#define UCHARFLOAT 3
#define CFLOAT 8
#define CUCHAR 5
#define CUCHARFLOAT 12

const int PNG_BYTES_TO_CHECK = 4;
const int ERROR = -1;

class ImageContent;
typedef ImageContent DARY;
typedef unsigned int uint;
typedef unsigned char uchar;

class ImageContent {
   
  private :
    uint x_size,y_size;
    uint tsize;
    void writePGM(const char *nom,unsigned char *buff, const char* comments);
    void write(const char *nom, const char* comments);
    void initFloat(uint ,uint);	   
    void initUChar(uint ,uint);	   
    void init3Float(uint ,uint);	   
    void init3UChar(uint ,uint);	   
    int buftype;
  
  public :
    float **fel;
    float **felr;
    float **felg;
    float **felb;
    char *filename;
    unsigned char **bel;
    unsigned char **belr;
    unsigned char **belg;
    unsigned char **belb;
    ImageContent(void){};   
    ImageContent(const char *);	
    ImageContent(const char *,int);
   
	  ImageContent(ImageContent *im);	   
	  ImageContent(uint y_size_in,uint x_size_in){initFloat( y_size_in, x_size_in);};	   
	  ImageContent(int y_size_in ,int x_size_in){initFloat((uint)y_size_in,(uint)x_size_in);};	   
	  ImageContent(uint ,uint, const char *);	   
	  ImageContent(uint ,uint, const char *, float);	   
	  ImageContent(uint y_size_in, uint x_size_in, float val){initFloat( y_size_in, x_size_in);set(val);};	   

	 ~ImageContent();

	  inline uint  const  x() const { return x_size;}
	  inline uint  const  y() const { return y_size;}
	  inline uint  const  size() const { return tsize;}
	  int const getType() const { return buftype;}
	  const char* name(){return filename;}
	  void write(const char *nom);
	  void writePNG(const char* name);
	  void writeR(const char *nom);
	  void writeG(const char *nom);
	  void writeB(const char *nom);
	  void RGB2xyY();
	  void RGB2lbrg();
	  void RGB2rgb();
	  void float2char();
	  void char2float();
	  void flipH();
	  void flipV();
	  void toGRAY();
	  void set(float);
	  void set(DARY*);
	  void set(const char *name){strcpy(filename,name);}
	  void normalize(float min_in, float max_in);
	  void normalize();
	  void scale(DARY *im_in, float scalex, float scaley);
          void scaleHalf(const DARY *im_in);
          float getValue(float x, float y);
	  void interpolate(DARY *sface, float m_x, float m_y, 
			   float scalex, float scaley, float angle);
	  void interpolate(DARY *im_in, float m_x, float m_y, float vec0x, float vec0y,
			   float vec1x, float vec1y);
          /*! copy from img so that coordinates img(x,y) maps to the middle of the patch */
	  void crop(DARY *img, int x, int y);
};
       


 
#endif
