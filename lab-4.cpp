#include <stdio.h>                                                                                                                          
#include <math.h>                                                                                                                           
#include <time.h>                                                                                                                           
#include <stdlib.h>                                                                                                                         
#include <mpi.h>                                                                                                                       
using namespace std;                                                                                                                      

int main(int argc, char* argv[]) {                                                                                                          

        int size, rank;                                                                                                                   
		
        MPI_Init(&argc, &argv);                                                                                                             
        MPI_Comm_size(MPI_COMM_WORLD, &size);                                                                                               
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                                                                               
		
        double *x = NULL, *y = NULL, *local_x, *local_y;                                                                                    
        double dot = 0, local_dot;                                                                                                          
        int n = 0, local_n;                                                                                                                 
        double start;                                                                                                                       
		
        if (rank == 0)                                                                                                                      
        {
			n = 1000000;                                                                                                                    
            x = (double *)calloc(n, sizeof(double));                                                                                        
            y = (double *)calloc(n, sizeof(double));                                                                                        
            for(int i = 0; i < n; i++) {                                                                                                    
                x[i] = y[i] = rand() / (RAND_MAX + 1.);                                                                                     
            }                                                                                                                               
            start = MPI_Wtime();                                                                                                            
        }                                                                                                                                   
		
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);                                                                                       
		
        local_n = n / size;                                                                                                                
		
        local_x = (double *)calloc(local_n, sizeof(double));                                                                                
        local_y = (double *)calloc(local_n, sizeof(double));                                                                                
		
        MPI_Scatter(x, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                                               
        MPI_Scatter(y, local_n, MPI_DOUBLE, local_y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                                             
		
        local_dot = 0;                                                                                                                      
        for (int i = 0; i < local_n; i++) {                                                                                                 
                local_dot += local_x[i] * local_y[i];                                                                                       
        }                                                                                                                                   
        free(local_x);                                                                                                                      
        free(local_y);                                                                                                                      
		
        if (rank == 0)
		{                                                                                                                                   
                for (int i = local_n * size; i < n; i++)                                                                                    
                        local_dot += x[i] * y[i];                                                                                           
        }                                                                                                                                  
		
        printf("Local dot of process #%d is: %f\n", rank, local_dot);                                                                       
        fflush(stdout);                                                                                                                     
		
        //double* dot_sum = (double *)malloc(sizeof(double) * ncpus);                                                                       
        //MPI_Allgather(&local_dot, 1, MPI_DOUBLE, dot_sum, 1, MPI_DOUBLE, MPI_COMM_WORLD);                                                 
        //then sum all elements of received dot_sums                                                                                        
		
        // recieve all 'local_dots' from processes and SUM them, then save it in root's 'dot'                                               
        MPI_Reduce(&local_dot, &dot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);                                                            
		
        if (rank == 0) {                                                                                                                    
                double end = MPI_Wtime();                                                                                                   
                double time = ((double)(end - start));                                                                                      
				
                printf("Processing time: %f, s\n", time);                                                                                   
                printf("Dot product: %f\n", dot);                                                                                           
        }                                                                                                                                  
		
        MPI_Finalize();                                                                                                                     
        return 0;                                                                                                                           
}