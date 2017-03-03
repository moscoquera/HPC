#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <mpi.h>
#include <math.h>


using namespace std;
const int max_val=10;

void generateArray(int* data, int size);
void solve(int* res, int* src1, int* src2, int size);


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
		solve(res,arr1,arr2,n);
	}else if (pid==0){
		int nnodes=np-1;
		int split=ceil((n*1.0)/nnodes);
		int* v1;
		int* v2;
		for(int i=0;i<nnodes-1;i++){
			v1=arr1+(split*i);
			v2=arr2+(split*i);
			int size=split*(i+1);
			MPI_Send(v1,size,MPI_INT,i+1,1,MPI_COMM_WORLD);
			MPI_Send(v2,size,MPI_INT,i+1,2,MPI_COMM_WORLD);
		}

		v1=arr1+(split*(nnodes-1));
		v2=arr2+(split*(nnodes-1));
		int size=n-(split*(nnodes-1));
		MPI_Send(v1,size,MPI_INT,nnodes,1,MPI_COMM_WORLD);
		MPI_Send(v2,size,MPI_INT,nnodes,2,MPI_COMM_WORLD);

		//datos procesados
		for(int nod=0;nod<nnodes;nod++){
			MPI_Status status;
			int number_amount;
			MPI_Probe(MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_INT, &number_amount);

			int* ptr;
			if (status.MPI_SOURCE==nnodes){ //es el ultimo nodo
				ptr=(int*)(res+(split*(nnodes-1)));
			}else{
				ptr=(int*)(res+(split*(status.MPI_SOURCE-1)));
			}
			//printf("%d %d %d %d %d\n",res,ptr,sizeof(unsigned long long int),split,(ptr-res));
			MPI_Recv(ptr,number_amount,MPI_INT,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}else if (pid>0){
		MPI_Status status;
		int number_amount;
		for(int tmp=0;tmp<2;tmp++){
			MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_INT, &number_amount);
			if(status.MPI_TAG==1){
				arr1 = (int*)malloc(sizeof(int)*number_amount);
				MPI_Recv(arr1,number_amount,MPI_INT,0,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}else if (status.MPI_TAG==2){
				arr2 = (int*)malloc(sizeof(int)*number_amount);
				MPI_Recv(arr2,number_amount,MPI_INT,0,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		res = (int*)malloc(sizeof(int)*number_amount);
		solve(res,arr1,arr2,number_amount);
		MPI_Send(res,number_amount,MPI_INT,0,3,MPI_COMM_WORLD);
	}

    
    if(pid==0){
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
	printf("\n");
    }
    MPI_Finalize();
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
