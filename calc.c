#include <stdlib.h>
#include <stdio.h>
#include<time.h>
#include<lapacke.h>
extern void LAPACK_dgeev( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );
extern void print_eigenvalues( char* desc, int n, double* wr, double* wi );
extern void print_eigenvectors( char* desc, int n, double* wi, double* v,
      int ldv );

double a[1500*1500];
double wr[1500*1500],wi[1500*1500],vl[1500*1500],vr[1500*1500];
int N,LDA,LDVL,LDVR;
int n, lda, ldvl, ldvr, info, lwork;
FILE *fptr;
FILE *fpt;
FILE *fp;
int main() {
	clock_t start, end;
    	double cpu_time_used;
    	start = clock();
    	fptr = fopen("matrix", "r");
    	fpt = fopen("valueseigen", "w");
    	fp = fopen("vectors","w");
	//N=getw(fptr);
	fscanf(fptr,"%d",&N);
	//scanf("%d",&N);
	LDA=N; LDVL=N; LDVR=N;
        n = N; lda = LDA; ldvl = LDVL; ldvr = LDVR;
        int i,j,k=0;
        double wkopt;
        double* work;
	//printf("Matrix is \n");
    	//fclose(fptr);
	for(int i=0;i<LDA*N;i++)
	{
	  double temp;
	  fscanf(fptr,"%lf",&temp);
		//a[i]=getw(fptr);
		a[i]=temp;
	}
        lwork = -1;
        LAPACK_dgeev( "Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Solve eigenproblem */
        LAPACK_dgeev( "Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         work, &lwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        print_eigenvalues( "Eigenvalues", n, wr, wi );
        /* Print left eigenvectors */
        print_eigenvectors( "Right eigenvectors", n, wi, vl, ldvl );
        /* Free workspace */
        free( (void*)work );
        fclose(fptr);
        fclose(fpt);
        fclose(fp);
	end = clock();
      	cpu_time_used = (double) (end - start)/(double)(CLOCKS_PER_SEC);
      	printf("Running time of code is %lf\n",cpu_time_used);
        exit( 0 );
} /* End of DGEEV Example */

/* Auxiliary routine: printing eigenvalues */
void print_eigenvalues( char* desc, int n, double* wr, double* wi ) {
        int j;
        //printf( "\n %s\n", desc );
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (double)0.0 ) {
         //printf( " %6.2f", wr[j] );
	 //putw(wr[j],fpt);
	 fprintf(fpt,"%lf ",wr[j]);	 
      } else {
         //printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
	 //putw(wr[j],fpt);
	 fprintf(fpt,"%lf ",wr[j]);
	 //putw(wi[j],fpt);
	 fprintf(fpt,"%lf ",wi[j]);
      }
   }
   //printf( "\n" );
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv ) {
        int i, j;
        //printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
            //printf( " %6.2f", v[i+j*ldv] );
            //putw(v[i+j*ldv],fp);
            fprintf(fp,"%lf ",v[i+j*ldv]);
            j++;
         } else {
            //printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            //putw(v[i+(j+1)*ldv],fp);
            fprintf(fp,"%lf ",v[i+j*ldv]);
            fprintf(fp,"%lf ",v[i+(j+1)*ldv]);
            //printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            //putw(-v[i+(j+1)*ldv],fp);
            fprintf(fp,"%lf ",v[i+j*ldv]);
            fprintf(fp,"%lf ",-v[i+(j+1)*ldv]);
            j += 2;
         }
      }
      //printf( "\n" );
   }
}
