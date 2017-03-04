#include "stdlib.h"
#include "stdio.h"
#include <math.h>
#include <cuda.h>

const int max_val=100;

void generateArray(float* data, int size);

__global__
void vectAddKernel(float* A, float* B, float* C, int n){
  int i = threadIdx.x+blockDim.x*blockIdx.x;
  if (i<n){
    *(C+i)=*(A+i)+*(B+i);
  }
}

void vectorAdd(float* A, float* B,float* C, int n){
  int size=sizeof(float)*n;
  

  float* d_A;
  float* d_B;
  float* d_C;
  
  int _sa = cudaMalloc((void**)(&d_A),size);
  int _sb = cudaMalloc((void**)(&d_B),size);
  int _sc = cudaMalloc((void**)(&d_C),size);

  int _cma=cudaMemcpy(d_A,A,size,cudaMemcpyHostToDevice);
  int _cmb=cudaMemcpy(d_B,B,size,cudaMemcpyHostToDevice);
  
  dim3 grid (ceil(n/256.0), 1, 1); 
  dim3 block (256, 1, 1);

  vectAddKernel<<<grid,block>>>(d_A,d_B,d_C,n);
  int _cmc=cudaMemcpy(C,d_C,size,cudaMemcpyDeviceToHost);
  

  cudaFree((void**)&d_A);
  cudaFree((void**)&d_B);
  cudaFree((void**)&d_C);

}


int main(int argc, char* argv[]){
    if (argc != 2){
      printf("Numero incorrecto de argumentos\n");
      return -1;
    }
    int n = atoi(argv[1]);
    
    float* arr1 = (float*)malloc(sizeof(float)*n);
    float* arr2 = (float*)malloc(sizeof(float)*n);
    float* res = (float*)malloc(sizeof(float)*n);
    generateArray(arr1,n);
    generateArray(arr2,n);

    vectorAdd(arr1,arr2,res,n);
    /*
    printf("Array 1:");
    for(int i=0;i<n;i++){
      printf(" %f",*(arr1+i));
    }
    printf("\n");
    
    printf("Array 2:");
    for(int i=0;i<n;i++){
      printf(" %f",*(arr2+i));
    }
    printf("\n");
    
    
    printf("Res:");
    for(int i=0;i<n;i++){
      printf(" %f",*(res+i));
    }
    printf("\n");*/

}

void generateArray(float* data, int size){
  for(int i=0;i<size;i++){
    *(data+i)=rand() % max_val;
  }
}

