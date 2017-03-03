#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <mpi.h>
#include <tgmath.h>


using namespace std;
const int max_val=10;

void generateArray(int* data, int size);
void solve(int* res, int* src1, int* src2, int from, int to);


int main(int argc, char* argv[]){
    if (argc != 2){
      printf("Numero incorrecto de argumentos\n");
      return -1;
    }
    int n = atoi(argv[1]);
	int pid=-1,np=-1;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&pid);
	
	MPI_Comm_size(MPI_COMM_WORLD,&np);    


    int* arr1;
    int* arr2;
    int* res;
	if (n<np){
		np=n;	
	}
	
	if (pid==0){
		arr1 = (int*)malloc(sizeof(int)*n);
    		arr2 = (int*)malloc(sizeof(int)*n);
		res = (int*)malloc(sizeof(int)*n);
		generateArray(arr1,n);
    		generateArray(arr2,n);		
	}
	if (np==1){
		solve(res,arr1,arr2,0,n);
	}else if (pid==0){
		int nnodes=np-1;
		int split=ceil((p*m*1.0)/nnodes);
		int* v1;
		int* v2;
		for(int i=0;i<nnodes-1;i++){
			v1=arr1+(split*i);
			v2=arr2+(split*i);
			int size=split*(i+1);
			MPI_Send(v1,size,MPI_INT,i+1,0,MPI_COMM_WORLD);
			MPI_Send(v2,size,MPI_INT,i+1,0,MPI_COMM_WORLD);
		}

		v1=arr1+(split*(nnodes-1));
		v2=arr2+(split*(nnodes-1));
		int size=n-(split*(nnodes-1));
		MPI_Send(v1,size,MPI_INT,i+1,0,MPI_COMM_WORLD);
		MPI_Send(v2,size,MPI_INT,i+1,0,MPI_COMM_WORLD);
	}else if (pid>0){

	}

    
    if(pid==0){
	/*printf("Array 1:");
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
    }
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
