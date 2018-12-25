#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
int arr[1500][1500];
int m[1500];
int M[1500][1500];
int bij[1500];
int I[1500][1500],Ad[1500][1500],Mt[1500][1500],B[1500][1500],C[1500][1500];
int n;
FILE *fptr;
FILE *fpt;
FILE *fp;
FILE *ft;
int main()
{
	clock_t start, end;
    	double cpu_time_used;
    	start = clock();
	fp=fopen("matrix","r");
	fptr=fopen("valueseigen","r");
	fpt=fopen("vectors","r");
	ft=fopen("ot_transpose","r");
	int i,j,k,dec;
	int n,s=0;
	fscanf(fp,"%d",&n);
	double values[n];
	double vectors[n][n];
	double matrix[n][n];
	double matmul[n][1];
	double lhs[n];
	double transpose[n][n];
	double finalvec[n][n];
	printf("Encrypted matrix was\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			fscanf(fp,"%lf",&matrix[i][j]);
			printf("%lf\t",matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	//printf("Eigen values of the matrix\n");
	for(i=0;i<n;i++)
	{
		fscanf(fptr, "%lf",&values[i]);
		//printf("%lf\t",values[i]);
	}
	
	printf("\n");
	//printf("Eigen Vectors are\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			fscanf(fpt,"%lf",&vectors[i][j]);
			//printf("%lf\t",vectors[i][j]);
		}
		//printf("\n");
	}
	//printf("\n");
	
	printf("Enter no of vectors you want to verify between 1 and %d\n",n);
	int choice,p;
	scanf("%d",&choice);
	for(p=1;p<=choice;p++)
	{
		int random=(int)rand()%n;
		double eva=values[random];
		double vec[n];
		for(j=0;j<n;j++)
		{
			vec[j]=vectors[j][random];
		}
		//printf("\n");
		for(i=0;i<n;i++)
		{
			for(j=0;j<1;j++)
			{
				matmul[i][j]=0;
			}
		}
		for(i=0;i<n;i++)
		{
			for(j=0;j<1;j++)
			{	
				double mtm=0.00;
				for(k=0;k<n;k++)
				{
					mtm+=(matrix[i][k]*vec[k]);
				}
				matmul[i][j]=mtm;
				//printf("%lf\t",matmul[i][j]);
				lhs[i]=matmul[i][j];
			}
		}
		//printf("\n");
		double t[n];
		int s=0;
		for(i=0;i<n;i++)
		{
			t[i]=eva*vec[i];
			//printf("%lf\t",t[i]);
		}
		//printf("\n");
		for(i=0;i<n;i++)
		{
			if(abs(t[i]-lhs[i])>0.1)
			{
				s=1;
			}
		}
	}
	fclose(fptr);
	fptr=fopen("valueseigen","r");
	if(s==0)
	{
		printf("Result verified! :)\n\n");
		printf("Original Matrix Eigenvalues\n");
		int alpha,s1;
		double z[n];
		FILE *fptrk;
		fptrk = fopen("keys", "r");
		fscanf(fptrk,"%d",&alpha);
		fscanf(fptrk,"%d",&s1);
		for(i=0;i<n;i++)
		{
			fscanf(fptr,"%lf",&z[i]);
			double ans=(z[i]-s1)/alpha;
			printf("%lf\t",ans);
		}
		printf("\n");
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				fscanf(ft,"%lf",&transpose[i][j]);
			}
		}
		printf("\nEigenVectors of Original Matrix\n");
		k=-1;
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				double s2=0.00;
				if(i==0)
				{
					if(transpose[i][n-1]==1.00){s2=s2+vectors[n-1][j];}
					else if(transpose[i][n-1]==-1.00){s2=s2-vectors[n-1][j];}
					else{printf("This is fuck\n");}
				}
				else
				{
					if(transpose[i][k]==1.00){s2=s2+vectors[k][j];}
					else if(transpose[i][k]==-1.00){s2=s2-vectors[k][j];}
					else{printf("This is fuck\n");}
				}
				finalvec[i][j]=s2;
				printf("%lf\t",finalvec[i][j]);
			}
			k=k+1;
			printf("\n");
		}
	}
	else
	{
		printf("Result obtained is wrong :(\n");
	}
	fclose(fp);
	fclose(fptr);
	fclose(fpt);
	fclose(ft);
	printf("\n");
	cpu_time_used = (double) (end - start)/(double)(CLOCKS_PER_SEC);
      	printf("Running time of code is %lf\n\n",cpu_time_used);
	return 0;
}
