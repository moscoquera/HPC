#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <cmath>

using namespace std;
using namespace cv;



int* filter(int* data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode);


int* uchartoint(uchar* data, int size){
	int* buff = (int*)malloc(sizeof(int)*size);
	for(int i=0;i<size;i++){
		*(buff+i)=(int)*(data+i);
	}
	return buff;
}


uchar* inttouchar(int* data, int size){
	uchar* buff = (uchar*)malloc(sizeof(uchar)*size);
	for(int i=0;i<size;i++){
		*(buff+i)=(unsigned char)*(data+i);
	}
	return buff;
}


int convolution(int* image,int channels, float *kernel, int kerneldim, int imagecols){
	float value=0;
	int pixel=0;
	int kernelmid = kerneldim%2==1?kerneldim/2:(kerneldim-1)/2;
	kernel=kernel+(((kerneldim+1)*kernelmid));
	for(int r = (-1*kerneldim/2); r<=(kerneldim/2);r++){
		for(int c = (-1*kerneldim/2); c<=(kerneldim/2);c++){
			pixel=*(image+(r*imagecols*channels)+(c*channels));
			value+=*(kernel+((r*-1)*kerneldim)+(c*-1))*pixel;
			
		}
	}
	return (int)value;

}


void normalize(int* data,int channels, int rows, int cols,int min, int max, int newMin, int newMax, int mode){
	int pixval=0;
	for(int r=1;r<rows-1;r++){
		for(int c=1;c<cols-1;c++){
			for(int cha=0;cha<channels;cha++){
				pixval=*(data+(r*cols*channels)+(c*channels)+cha);
				if (mode==1){
					*(data+(r*cols*channels)+(c*channels)+cha)=(pixval-min)*((newMax-newMin*1.0)/(max-min))+newMin;
				}else{
					*(data+(r*cols*channels)+(c*channels)+cha)=pixval>newMax?newMax:pixval<newMin?newMin:pixval;
				}
	    		}
		}
	  }
}


uchar * edge1(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,0,-1},
		{0,0,0},
		{-1,0,1}
	};
	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,0,0);
	return inttouchar(bf,rows*cols*channels);
}


uchar * edge2(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{0,1,0},
		{1,-4,1},
		{0,1,0}
	};

	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,0,0);
	return inttouchar(bf,rows*cols*channels);
}


uchar * edge3(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{-1,-1,-1},
		{-1,8,-1},
		{-1,-1,-1},
	};
	
	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,0,0);
	return inttouchar(bf,rows*cols*channels);
}

uchar * sharpen(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{0,-1,0},
		{-1,5,-1},
		{0,-1,0},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,0,0);
	return inttouchar(bf,rows*cols*channels);
}

uchar * boxblur(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,1,1},
		{1,1,1},
		{1,1,1},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,1,0);
	return inttouchar(bf,rows*cols*channels);
}

uchar * gaussianblur(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,2,1},
		{2,4,2},
		{1,2,1},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	bf=filter(bf,channels,rows,cols,(float*)kernel,3,1,0);
	return inttouchar(bf,rows*cols*channels);
}

int * _sobelx(uchar* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{1,0,-1},
		{2,0,-2},
		{1,0,-1},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	return filter(bf,channels,rows,cols,(float*)kernel,3,0,mode);
	
}
uchar * sobelx(uchar* data,int channels, int rows,int cols){
	return inttouchar(_sobelx(data,channels,rows,cols,0),rows*cols*channels);
}

int * _sobely(uchar* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{1,2,1},
		{0,0,0},
		{-1,-2,-1},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	return filter(bf,channels,rows,cols,(float*)kernel,3,0,mode);
	
}

uchar * sobely(uchar* data,int channels, int rows,int cols){
	return inttouchar(_sobely(data,channels,rows,cols,0),rows*cols*channels);
}


uchar * sobel(uchar* data,int channels, int rows,int cols){
	 int * filterx =  _sobelx(data,channels,rows,cols,-1);
	 int * filtery =  _sobely(data,channels,rows,cols,-1);
	 int * output = (int*)malloc(sizeof(int)*rows*cols*channels);
	int val;
	int min=INT_MAX;
	int max=INT_MIN;
	for(int i=0;i<rows*cols*channels;i++){
		val=sqrt((*(filterx+i))*(*(filterx+i))+(*(filtery+i))*(*(filtery+i)));
		*(output+i)=val;
		if (val<min) min=val;
		if (val>max) max=val;
	}
	normalize(output,channels,rows,cols,min,max,0,255,1);	
	return inttouchar(output,rows*cols*channels);
	 
}


int * _sobelx10(uchar* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{3,0,-3},
		{10,0,-10},
		{3,0,-3},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	return filter(bf,channels,rows,cols,(float*)kernel,3,0,mode);
}

int * _sobely10(uchar* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{3,10,3},
		{0,0,0},
		{-3,-10,-3},
	};

	int* bf = uchartoint(data,rows*cols*channels);
	return filter(bf,channels,rows,cols,(float*)kernel,3,0,mode);
}

uchar * sobel10(uchar* data,int channels, int rows,int cols){
	  int * filterx =  _sobelx10(data,channels,rows,cols,-1);
	 int * filtery =  _sobely10(data,channels,rows,cols,-1);
	 int * output = (int*)malloc(sizeof(int)*rows*cols*channels);
	int val;
	int min=INT_MAX;
	int max=INT_MIN;
	for(int i=0;i<rows*cols*channels;i++){
		val=sqrt((*(filterx+i))*(*(filterx+i))+(*(filtery+i))*(*(filtery+i)));
		*(output+i)=val;
		if (val<min) min=val;
		if (val>max) max=val;
	}
	normalize(output,channels,rows,cols,min,max,0,255,1);	
	return inttouchar(output,rows*cols*channels);
	 
}

uchar * sobely10(uchar* data,int channels, int rows,int cols){
	return inttouchar(_sobely10(data,channels,rows,cols,0),rows*cols*channels);
}

uchar * sobelx10(uchar* data,int channels, int rows,int cols){
	return inttouchar(_sobelx10(data,channels,rows,cols,0),rows*cols*channels);
}


int* filter(int* data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode){
	int* buff;
	buff=(int*)malloc(sizeof(int)*channels*rows*cols);
	for(int i=0;i<cols*rows*channels;i++){
		*(buff+i)=0; //clean up
	}

	int max=INT_MIN;
	int min=INT_MAX;
	int newMax=INT_MIN;
	int newMin=INT_MAX;

	if (kernelNormalize){
		float sumKernel[2]={0,0};
		float kernelVal;
		for(int i=0;i<9;i++){
			kernelVal=*((float*)kernel+i);
			sumKernel[kernelVal>=0]+=kernelVal;
		}
		if (kernelNormalize==1){
			sumKernel[0]+=sumKernel[1];
			sumKernel[1]=sumKernel[0];
		}
		for(int i=0;i<9;i++){
			kernelVal=*((float*)kernel+i);
			*((float*)kernel+i)=kernelVal/sumKernel[kernelVal>=0];
		}
	}
	int value,pixel;
	  for(int r=1;r<rows-1;r++){
	      for(int c=1;c<cols-1;c++){
		for(int cha=0;cha<channels;cha++){
			value=convolution((data+(r*cols*channels)+(c*channels)+cha),channels,kernel,3,cols);
			if (value>max){
				max=value;
			}

			if (value<min){
				min=value;	
			}
			*(buff+(r*cols*channels)+(c*channels)+cha)=value;
			
			pixel=*(data+(r*cols*channels)+(c*channels)+cha);
			if (pixel<newMin){
				newMin=pixel;			
			}
			if (pixel>newMax){
				newMax=pixel;			
			}
			
		    }
		}
	  }
	printf("%d %d %d %d\n",min,max,newMin,newMax);
	if (outputNormalizationMode>=0){
		  normalize(buff,channels,rows,cols,min,max,newMin,newMax,outputNormalizationMode);
	}
	  return buff;

}




int main(int argc, char** argv){


	if (argc<3){
		cout<<"./nombre imagen filtro"<<endl;
		return 0;
	}


	char* nfiltro=*(argv+2);
	uchar* (*filtro)(uchar*,int,int,int)=0;

	if(strcmp(nfiltro,"sobel")==0) filtro=sobel;
	if(strcmp(nfiltro,"sobelx")==0) filtro=sobelx;
	if(strcmp(nfiltro,"sobely")==0) filtro=sobely;
	if(strcmp(nfiltro,"sobel10")==0) filtro=sobel10;
	if(strcmp(nfiltro,"sobelx10")==0) filtro=sobelx10;
	if(strcmp(nfiltro,"sobely10")==0) filtro=sobely10;
	if(strcmp(nfiltro,"edge1")==0) filtro=edge1;
	if(strcmp(nfiltro,"edge2")==0) filtro=edge2;
	if(strcmp(nfiltro,"edge3")==0) filtro=edge3;
	if(strcmp(nfiltro,"boxblur")==0) filtro=boxblur;
	if(strcmp(nfiltro,"gaussianblur")==0) filtro=gaussianblur;
	if(strcmp(nfiltro,"sharpen")==0) filtro=sharpen;
	if (filtro==0){
		cout<<"metodo erroneo"<<endl;
		return 1;
	}

	Mat image;
	image = imread(*(argv+1), CV_LOAD_IMAGE_COLOR);
	Mat m1;

    
	if(! image.data )                              // Check for invalid input
	{
		cout <<  "Could not open or find the image" << std::endl ;
		return -1;
	}


	m1 = Mat (image);
	m1.data=filtro(image.data,3,image.rows,image.cols);
	namedWindow( "original", WINDOW_AUTOSIZE );
	imshow( "original", image );             
	namedWindow( "filter", WINDOW_AUTOSIZE );
	imshow( "filter", m1 );             

    waitKey();                                        // Wait for a keystroke in the window
  return 0;
}


