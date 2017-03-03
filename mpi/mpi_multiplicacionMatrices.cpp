#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <mpi.h>
#include <tgmath.h>


//using namespace std;
const int max_val=100;

void generateArray(int* data, int size);
void solve(unsigned long long int* res, int* src1, int* src2, int n,int p, int from,int to);


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
	
    if (argc < 4){
      printf("Numero incorrecto de argumentos\n");
      return -1;
    }
	int m = atoi(argv[1]);
	int n = atoi(argv[2]);
	int p = atoi(argv[3]);

	int pid=-1,np=-1;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&pid);
	
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	int* mat1;
	int* mat2;
	unsigned long long int* res;

	if (m*p<np){
		np=m*p;	
	}
	
	if (pid==0){
		mat1 = (int*)malloc(sizeof(int)*m*n);
		mat2 = (int*)malloc(sizeof(int)*p*n);
		res = (unsigned long long int*)malloc(sizeof(unsigned long long int)*m*p);
		generateArray(mat1,n*m);
		generateArray(mat2,n*p);
		
	}
	if (np==1){
		solve(res,mat1,mat2,n,p,0,m*p);
	}else if (pid==0){
		//envia mat1 y mat2
		for(int i=1;i<np;i++){
			MPI_Send(mat1,m*n,MPI_INT,i,1,MPI_COMM_WORLD);
			MPI_Send(mat2,p*n,MPI_INT,i,2,MPI_COMM_WORLD);
		}
		int nnodes=np-1;
		int split=ceil((p*m*1.0)/nnodes);
		int vfrom;
		int vto;
		for(int i=0;i<nnodes-1;i++){
			vfrom=split*i;
			vto=split*(i+1);
			MPI_Send(&vfrom,1,MPI_INT,i+1,3,MPI_COMM_WORLD);
			MPI_Send(&vto,1,MPI_INT,i+1,3,MPI_COMM_WORLD);
		}
		vfrom=split*(nnodes-1);
		vto=m*p;
		MPI_Send(&vfrom,1,MPI_INT,nnodes,3,MPI_COMM_WORLD);
		MPI_Send(&vto,1,MPI_INT,nnodes,3,MPI_COMM_WORLD);

		//datos procesados
		
		for(int nod=0;nod<nnodes;nod++){
			MPI_Status status;
			int number_amount;
			MPI_Probe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_UNSIGNED_LONG_LONG, &number_amount);

			unsigned long long int* ptr;
			if (status.MPI_SOURCE==nnodes){ //es el ultimo nodo
				ptr=(unsigned long long int*)(res+(split*(nnodes-1)));
			}else{
				ptr=(unsigned long long int*)(res+(split*(status.MPI_SOURCE-1)));
			}
			//printf("%d %d %d %d %d\n",res,ptr,sizeof(unsigned long long int),split,(ptr-res));
			MPI_Recv(ptr,number_amount,MPI_UNSIGNED_LONG_LONG,MPI_ANY_SOURCE,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			
		}
		

	}else if (pid>0){
			MPI_Status status;
			int number_amount;
			for(int tmp=0;tmp<2;tmp++){
				MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status, MPI_INT, &number_amount);
				if(status.MPI_TAG==1){
					mat1 = (int*)malloc(sizeof(int)*number_amount);
					MPI_Recv(mat1,number_amount,MPI_INT,0,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}else if (status.MPI_TAG==2){
					mat2 = (int*)malloc(sizeof(int)*number_amount);
					MPI_Recv(mat2,number_amount,MPI_INT,0,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				
			}
			
			//recibe que datos procesar
			int from,to;
			MPI_Recv(&from,1,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&to,1,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//printf("node: %d from: %d to: %d\n",pid,from,to);
			res = (unsigned long long int*)malloc(sizeof(unsigned long long int)*(to-from));
			solve(res,mat1,mat2,n,p,from,to);
			MPI_Send(res,(to-from),MPI_UNSIGNED_LONG_LONG,0,4,MPI_COMM_WORLD);

		}

    /*
    int* mat1 = (int*)malloc(sizeof(int)*m*n);
    int* mat2 = (int*)malloc(sizeof(int)*p*n);
    unsigned long long int* res = (unsigned long long int*)malloc(sizeof(unsigned long long int)*m*p);
    generateArray(mat1,n*m);
    generateArray(mat2,n*p);
    solve(res,mat1,mat2,m,n,p);
    */

    if(pid==0){
    	imprimir((char*)"Mat 1:",m,n,mat1);
	imprimir((char*)"Mat 2:",n,p,mat2);
	imprimirllu((char*)"Res :",m,p,res);
    }
    
    MPI_Finalize();
    
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

void solve(unsigned long long int* res, int* src1, int* src2, int n,int p, int from,int to){
	for (int index=from;index<to;index++){
		int c=index%p;
		int f=index/p;
		*(res+(index-from))=generateCell(src1,src2,n,p,f,c);		
	}

}
