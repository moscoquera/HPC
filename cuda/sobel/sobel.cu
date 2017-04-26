#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <cuda.h>


using namespace std;
using namespace cv;



int* _filter(int* data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode);
uchar * filter(uchar * data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode);
__device__ int getGlobalIdx_3D_3D();
__device__ int getblockthreadIdx();


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


__global__
void convolution(int* data,int* buff,float* kernel,int* outputvars,int rows,int cols,int channels,int kerneldim){
	int idx = getGlobalIdx_3D_3D();
	int kernelmid;
	extern __shared__ float sharedKernel[];
	float *kernelCenter;

	if (getblockthreadIdx()<kerneldim*kerneldim){
		*(sharedKernel+getblockthreadIdx())=*(kernel+getblockthreadIdx());
	}

	__syncthreads();
/*	
	if (getblockthreadIdx()<kerneldim*kerneldim){
		printf("%d %f\n",getblockthreadIdx(),*(sharedKernel+getblockthreadIdx()));
	}

	__syncthreads();
*/
	kernelmid = kerneldim%2==1?kerneldim/2:(kerneldim-1)/2;
	kernelCenter=sharedKernel+(((kerneldim+1)*kernelmid));
	int row = idx / (cols*channels);
	int col = (idx%(cols*channels))/channels;
	float value=0;
	int pixel=0;
	float kernelVal=0;
	int pixelmin=INT_MAX,pixelmax=INT_MIN;
	int kernelmidHalf=(kerneldim/2);
	if (col>0 && row>0 && row<rows-1 && col<cols-1){
		data = data+idx;
		//r<=(kernelmidHalf) no funciona, no sé porque, pero cuda y yo tenemos un problema.
		for(int r = (-1*kernelmidHalf); r<(kernelmidHalf+1);r++){
			for(int c = -1*kernelmidHalf; c<(kernelmidHalf+1);c++){
				pixel=*(data+(r*cols*channels)+(c*channels));
				kernelVal=*(kernelCenter+(r*-1*kerneldim)+(c*-1));
				value+=kernelVal*pixel;
				if (pixel<pixelmin){
					pixelmin=pixel;			
				}
				if (pixel>pixelmax){
					pixelmax=pixel;			
				}
				
				
			}

		}
		*(buff+idx)=value;

		atomicMin(outputvars,value);	
		atomicMax(outputvars+1,value);
		atomicMin(outputvars+2,pixelmin);	
		atomicMax(outputvars+3,pixelmax);
	
	}	
	//__syncthreads();

	/*if (col>0 && row>0 && row<rows-1 && col<cols-1 && getblockthreadIdx()==0){
		printf("%d %d %d %d\n",*(outputvars),*(outputvars+1),*(outputvars+2),*(outputvars+3));
	}*/

}

__global__
void normalize(int* data,int channels, int rows, int cols,int min, int max, int newMin, int newMax, int mode){
	int pixval=0;
	int i = getGlobalIdx_3D_3D();
	int row = i / (cols*channels);
	int col = (i%(cols*channels))/channels;
	if (row>0 && col>0 && row<rows-1 && col<cols-1){
		pixval=*(data+i);
		if (mode==1){
			*(data+i)=(pixval-min)*((newMax-newMin*1.0)/(max-min))+newMin;
		}else{
			*(data+i)=pixval>newMax?newMax:pixval<newMin?newMin:pixval;
		}
	}
//	__syncthreads();

}


uchar * edge1(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,0,-1},
		{0,0,0},
		{-1,0,1}
	};

	return filter(data,channels,rows,cols,(float*)kernel,3,0,0);
}


uchar * edge2(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{0,1,0},
		{1,-4,1},
		{0,1,0}
	};
	
 	return filter(data,channels,rows,cols,(float*)kernel,3,0,0);

}


uchar * edge3(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{-1,-1,-1},
		{-1,8,-1},
		{-1,-1,-1},
	};
	
	return filter(data,channels,rows,cols,(float*)kernel,3,0,0);
}

uchar * sharpen(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{0,-1,0},
		{-1,5,-1},
		{0,-1,0},
	};

	return filter(data,channels,rows,cols,(float*)kernel,3,0,0);

}

uchar * boxblur(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,1,1},
		{1,1,1},
		{1,1,1},
	};

	return filter(data,channels,rows,cols,(float*)kernel,3,1,0);
	
}

uchar * gaussianblur(uchar* data,int channels, int rows,int cols){
	float kernel[3][3]={
		{1,2,1},
		{2,4,2},
		{1,2,1},
	};

	return filter(data,channels,rows,cols,(float*)kernel,3,1,0);

}

int * _sobelx(int* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{1,0,-1},
		{2,0,-2},
		{1,0,-1},
	};
	float * d_kernel;
	cudaMalloc(&d_kernel,sizeof(float)*3*3);
	cudaMemcpy(d_kernel,(float*)kernel,sizeof(float)*3*3,cudaMemcpyHostToDevice);
	int* res = _filter(data,channels,rows,cols,d_kernel,3,0,mode);
	cudaFree(d_kernel);
	return res;
	
	
	
}
uchar * sobelx(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data;
	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);
	int* d_output = _sobelx(d_data,channels,rows,cols,0);
	int* output = (int*)malloc(sizeof(int)*rows*cols*channels);
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(d_data);
	cudaFree(d_output);
	free(output);
	return out;
}

int * _sobely(int* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{1,2,1},
		{0,0,0},
		{-1,-2,-1},
	};
	float * d_kernel;
	cudaMalloc(&d_kernel,sizeof(float)*3*3);
	cudaMemcpy(d_kernel,(float*)kernel,sizeof(float)*3*3,cudaMemcpyHostToDevice);
	int* res = _filter(data,channels,rows,cols,d_kernel,3,0,mode);
	cudaFree(d_kernel);
	return res;
	
}

int * _sobelx10(int* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{3,0,-3},
		{10,0,-10},
		{3,0,-3},
	};

	float * d_kernel;
	cudaMalloc(&d_kernel,sizeof(float)*3*3);
	cudaMemcpy(d_kernel,(float*)kernel,sizeof(float)*3*3,cudaMemcpyHostToDevice);
	int* res = _filter(data,channels,rows,cols,d_kernel,3,0,mode);
	cudaFree(d_kernel);
	return res;
	
}

int * _sobely10(int* data,int channels, int rows,int cols, int mode){
	float kernel[3][3]={
		{3,10,3},
		{0,0,0},
		{-3,-10,-3},
	};

	float * d_kernel;
	cudaMalloc(&d_kernel,sizeof(float)*3*3);
	cudaMemcpy(d_kernel,(float*)kernel,sizeof(float)*3*3,cudaMemcpyHostToDevice);
	int* res = _filter(data,channels,rows,cols,d_kernel,3,0,mode);
	cudaFree(d_kernel);
	return res;
	
}


uchar * sobely(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data;
	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);
	int* d_output = _sobely(d_data,channels,rows,cols,0);
	int* output = (int*)malloc(sizeof(int)*rows*cols*channels);
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(d_data);
	cudaFree(d_output);
	free(output);
	return out;
}

__global__
void sobelKernel(int *a, int*b,int* output,int* outputvars,int n){
	int i = getGlobalIdx_3D_3D();
	
	if (i>=n){return;}
	int val=sqrtf((*(a+i))*(*(a+i))+(*(b+i))*(*(b+i)));
	*(output+i)=val;
	atomicMin(outputvars,val);	
	atomicMax(outputvars+1,val);
}

uchar * sobel(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data,*minmaxs;
	int * d_output,*output;

	cudaMalloc(&minmaxs,sizeof(int)*2);

	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);


	cudaMalloc(&d_output,sizeof(int)*rows*cols*channels);
	output = (int*)malloc(sizeof(int)*rows*cols*channels);
	

	int * filterx =  _sobelx(d_data,channels,rows,cols,-1);
	int * filtery =  _sobely(d_data,channels,rows,cols,-1);
	cudaMemset(minmaxs,INT_MAX,1);
	cudaMemset(minmaxs+1,INT_MIN,1);
	sobelKernel<<<ceil((rows*cols*channels)/256.0),256>>>(filterx,filtery,d_output,minmaxs,rows*cols*channels);
	int* tmpMinMax = (int*)malloc(sizeof(int)*2);	
	cudaMemcpy(tmpMinMax,minmaxs,sizeof(int)*2, cudaMemcpyDeviceToHost);

	normalize<<<ceil((rows*cols*channels)/256.0),256>>>(d_output,channels,rows,cols,*(tmpMinMax),*(tmpMinMax+1),0,255,1);	
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	//printf("%d %d %d %d\n",*(tmpMinMax),*(tmpMinMax+1),0,255);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(minmaxs);
	cudaFree(d_data);
	cudaFree(d_output);
	free(datai);
	free(output);
	cudaFree(filterx);
	cudaFree(filtery);
	free(tmpMinMax);
	return out;


}


uchar * sobel10(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data,*minmaxs;
	int * d_output,*output;

	cudaMalloc(&minmaxs,sizeof(int)*2);

	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);


	cudaMalloc(&d_output,sizeof(int)*rows*cols*channels);
	output = (int*)malloc(sizeof(int)*rows*cols*channels);
	

	int * filterx =  _sobelx10(d_data,channels,rows,cols,-1);
	int * filtery =  _sobely10(d_data,channels,rows,cols,-1);
	cudaMemset(minmaxs,INT_MAX,1);
	cudaMemset(minmaxs+1,INT_MIN,1);
	sobelKernel<<<ceil((rows*cols*channels)/256.0),256>>>(filterx,filtery,d_output,minmaxs,rows*cols*channels);
	int* tmpMinMax = (int*)malloc(sizeof(int)*2);	
	cudaMemcpy(tmpMinMax,minmaxs,sizeof(int)*2, cudaMemcpyDeviceToHost);

	normalize<<<ceil((rows*cols*channels)/256.0),256>>>(d_output,channels,rows,cols,*(tmpMinMax),*(tmpMinMax+1),0,255,1);	
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	//printf("%d %d %d %d\n",*(tmpMinMax),*(tmpMinMax+1),0,255);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(minmaxs);
	cudaFree(d_data);
	cudaFree(d_output);
	free(datai);
	free(output);
	cudaFree(filterx);
	cudaFree(filtery);
	free(tmpMinMax);
	return out;
	 
}

uchar * sobely10(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data;
	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);
	int* d_output = _sobely10(d_data,channels,rows,cols,0);
	int* output = (int*)malloc(sizeof(int)*rows*cols*channels);
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(d_data);
	cudaFree(d_output);
	free(output);
	return out;
}

uchar * sobelx10(uchar* data,int channels, int rows,int cols){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data;
	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);
	int* d_output = _sobelx10(d_data,channels,rows,cols,0);
	int* output = (int*)malloc(sizeof(int)*rows*cols*channels);
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	uchar* out = inttouchar(output,rows*cols*channels);
	cudaFree(d_data);
	cudaFree(d_output);
	free(output);
	return out;
	
}


__global__
void kernelNormAdd(float* kernel,float* output, int kernelNormalize){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	float kernelVal=*((float*)kernel+i);

	atomicAdd(output+(kernelVal>=0),kernelVal);
	__syncthreads();
	if (kernelNormalize==1){
		*(kernel+i)=kernelVal/(*output+*(output+1));
	}else{
		*(kernel+i)=kernelVal/(*(output+(kernelVal>=0)));
	}
	__syncthreads();
	
}


int* _filter(int* data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode){
	int* buff,*minmaxs;
	cudaMalloc(&buff,sizeof(int)*channels*rows*cols);
	cudaMalloc(&minmaxs,sizeof(int)*4);
	cudaMemset(buff,0,sizeof(int)*channels*rows*cols);
	if (kernelNormalize){
		float* sumKernel;
		cudaMalloc(&sumKernel,sizeof(float)*2);
		cudaMemset(sumKernel,0,sizeof(float)*2);
		kernelNormAdd<<<1,9>>>(kernel,sumKernel,kernelNormalize);
		cudaFree(sumKernel);
	}
	int N = rows*cols*channels;
	int ssize = (sizeof(float)*kerneldim*kerneldim);
	cudaMemset(minmaxs,INT_MAX,1);
	cudaMemset(minmaxs+1,INT_MIN,1);
	cudaMemset(minmaxs+2,INT_MAX,1);
	cudaMemset(minmaxs+3,INT_MIN,1);
	printf("%f\n",ceil(N/512.0));
	convolution<<<ceil(N/512.0),512,ssize>>>(data,buff,kernel,minmaxs,rows,cols,channels,kerneldim);
	cudaError_t err=cudaGetLastError();
	if ( cudaSuccess !=  err ){
	    printf( "Error!\n" );
            printf("GPUassert: %s\n", cudaGetErrorString(err));
	}
	if (outputNormalizationMode>=0){
		  int* tmpMinMax = (int*)malloc(sizeof(int)*4);	
		  cudaMemcpy(tmpMinMax,minmaxs,sizeof(int)*4, cudaMemcpyDeviceToHost);
		  //printf("%d %d %d %d\n",*(tmpMinMax),*(tmpMinMax+1),*(tmpMinMax+2),*(tmpMinMax+3));
		  normalize<<<ceil(N/256),256>>>(buff,channels,rows,cols,*(tmpMinMax),*(tmpMinMax+1),*(tmpMinMax+2),*(tmpMinMax+3),outputNormalizationMode);
		free(tmpMinMax);
	}
	cudaFree(minmaxs);
	return buff;
}

uchar* filter(uchar* data,int channels, int rows,int cols,float *kernel,int kerneldim, int kernelNormalize, int outputNormalizationMode){
	int* datai = uchartoint(data,channels*rows*cols);
	int * d_data;
	cudaMalloc(&d_data,sizeof(int)*channels*rows*cols);
	cudaMemcpy(d_data,datai,sizeof(int)*channels*rows*cols,cudaMemcpyHostToDevice);
	float * d_kernel;
	cudaMalloc(&d_kernel,sizeof(float)*3*3);
	cudaMemcpy(d_kernel,kernel,sizeof(float)*3*3,cudaMemcpyHostToDevice);
	int* d_output = _filter(d_data,channels,rows,cols,d_kernel,kerneldim,kernelNormalize,outputNormalizationMode);
	int* output = (int*)malloc(sizeof(int)*rows*cols*channels);
	cudaMemcpy(output,d_output,sizeof(int)*rows*cols*channels, cudaMemcpyDeviceToHost);
	uchar* out = inttouchar(output,rows*cols*channels);

	cudaFree(d_data);
	cudaFree(d_output);
	free(datai);
	free(output);
	return out;
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


__device__ int getGlobalIdx_3D_3D()
{
	int blockId = blockIdx.x 
			 + blockIdx.y * gridDim.x 
			 + gridDim.x * gridDim.y * blockIdx.z; 
	int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
			  + (threadIdx.z * (blockDim.x * blockDim.y))
			  + (threadIdx.y * blockDim.x)
			  + threadIdx.x;
	return threadId;
}


__device__ int getblockthreadIdx(){
	return (threadIdx.z * (blockDim.x * blockDim.y))
			  + (threadIdx.y * blockDim.x)
			  + threadIdx.x;
}



