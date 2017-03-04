#include "stdlib.h"
#include "stdio.h"
#include "string.h"

using namespace std;
const int max_val=100;

void generateArray(int* data, int size);
void solve(int* res, int* src1, int* src2, int size);


int main(int argc, char* argv[]){
    if (argc != 2){
      printf("Numero incorrecto de argumentos\n");
      return -1;
    }
    int n = atoi(argv[1]);
    
    int* arr1 = (int*)malloc(sizeof(int)*n);
    int* arr2 = (int*)malloc(sizeof(int)*n);
    int* res = (int*)malloc(sizeof(int)*n);
    generateArray(arr1,n);
    generateArray(arr2,n);
    solve(res,arr1,arr2,n);
    /*
    printf("Array 1:");
    for(int i=0;i<n;i++){
      printf(" %d",*(arr1+i));
    }
    printf("\n");
    
    printf("Array 2:");
    for(int i=0;i<n;i++){
      printf(" %d",*(arr2+i));
    }
    printf("\n");
    
    
    printf("Res:");
    for(int i=0;i<n;i++){
      printf(" %d",*(res+i));
    }
    printf("\n");*/
    
    return 0;
}

void generateArray(int* data, int size){
  for(int i=0;i<size;i++){
    *(data+i)=rand() % max_val;
  }
}

void solve(int* res, int* src1, int* src2, int size){
  for(int i=0;i<size;i++){
    *(res+i)=(*(src1+i)+*(src2+i));
  }
}
