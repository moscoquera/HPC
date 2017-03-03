#include "stdlib.h"
#include "stdio.h"
#include "string.h"

using namespace std;
const int max_val=100;

void generateArray(int* data, int size);
void solve(unsigned long long int* res, int* src1, int* src2, int m, int n,int p);


void imprimir(char* name, int rows, int cols, int* data){

    printf(name);
    printf("\n");
    for(int f=0;f<rows;f++){
      for(int c=0;c<cols;c++){
      printf(" %d",*(data+(f*cols)+c));
    }
      printf("\n");
    }
    printf("\n");

}

void imprimirllu(char* name, int rows, int cols, unsigned long long int* data){

    printf(name);
    printf("\n");
    for(int f=0;f<rows;f++){
      for(int c=0;c<cols;c++){
      printf(" %llu",*(data+(f*cols)+c));
    }
      printf("\n");
    }
    printf("\n");

}


int main(int argc, char* argv[]){
    if (argc != 4){
      printf("Numero incorrecto de argumentos\n");
      return -1;
    }
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int p = atoi(argv[3]);
    
    int* mat1 = (int*)malloc(sizeof(int)*m*n);
    int* mat2 = (int*)malloc(sizeof(int)*p*n);
    unsigned long long int* res = (unsigned long long int*)malloc(sizeof(unsigned long long int)*m*p);
    generateArray(mat1,n*m);
    generateArray(mat2,n*p);
    solve(res,mat1,mat2,m,n,p);
    
    imprimir((char*)"Mat 1:",m,n,mat1);
    imprimir((char*)"Mat 2:",n,p,mat2);
    imprimirllu((char*)"Res :",m,p,res);
    
    
    
    
    return 0;
}

void generateArray(int* data, int size){
  for(int i=0;i<size;i++){
    *(data+i)= rand() % max_val;
  }
}

int generateCell(int* src1, int* src2, int n,int p, int i,int j){
	int res=0;
	for(int k=0;k<n;k++){
		res+=*(src1+(i*n)+k)*(*(src2+(k*p)+j));
	}
	return res;
}

void solve(unsigned long long int* res, int* src1, int* src2, int m, int n,int p){
    for(int f=0;f<m;f++){
      	for(int c=0;c<p;c++){
      		*(res+(f*p)+c)=generateCell(src1,src2,n,p,f,c);
    	}
   }
}
