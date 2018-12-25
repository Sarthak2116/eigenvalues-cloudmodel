#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
double arr[1500][1500];
double m[1500];
double M[1500][1500];
double bij[1500];
double I[1500][1500],Ad[1500][1500],Mt[1500][1500],B[1500][1500],C[1500][1500];
int n;
FILE *fptr;
FILE *fpt;
FILE *fp;
int input()
{
	int i,j,k,dec;
	printf("Enter the order of matrix\n");
	scanf("%d",&n);
	//putw(n, fptr);
	fprintf(fptr,"%d ",n);
	printf("Do you want to enter the elements manually or want to automate this process?\n1 for Manual Entering Elements\n2 for Automating the Process\n");
	scanf("%d",&dec);
	srand(time(0));
	if(dec==1)
	{
		printf("Enter the elements of the matrix\n");
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				scanf("%lf",&arr[i][j]);
			}
		}
		return 1;
	}
	else if(dec==2)
	{
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				arr[i][j]=(rand()%1000)+1;
			}
		}
		return 1;
	}
	else
	{
		printf("Wrong choice\n");
		return 0;
	}
}
void keygen()
{
	int i,j;
	double range[]={-1.0,1.0};
	srand(time(0));
	for(i=0;i<n;i++)
	{
		if(i!=(n-1))
		{
			bij[i]=i+1;
		}
		else
		{
			bij[i]=0;
		}
	}
	for(i=0;i<n;i++)
	{
		int u=rand()%2;
		m[i]=range[u];
		
		for(j=0;j<n;j++)
		{
			if(bij[i]==j)
			{
				M[i][j]=m[i];
			}
			else
			{
				M[i][j]=0.00;
			}
		}
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			printf("%lf\t",M[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
void encrypt()
{
	int i,j,k;
	srand(time(0));
	int upper=3;
	int  lower=1;
	int alpha = (rand() % (upper - lower + 1)) + lower; 
	int s = (rand() % (upper - lower + 1)) + lower; 
	printf("alpha=%d\n",alpha);
	printf("s=%d\n",s);
	//putw(alpha, fpt);
	fprintf(fpt,"%d ",alpha);
	//putw(s, fpt);
	fprintf(fpt,"%d ",s);
	for(i=0;i<n;i++)
	{
	I[i][i]=1;
	}
	//A'=alpha(A)+s(I)
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			Ad[i][j] = ( (alpha*arr[i][j])+ (s*I[i][j]) );
		}
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			Mt[i][j]=M[j][i];
		}
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			fprintf(fp,"%lf ",Mt[i][j]);
		}
	}
	//Transpose(n);
	k=1;	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			double s1=0.00;
			if(i==n-1)
			{
				if(M[0][i]==1.00){s1=s1+Ad[0][j];}
				else if(M[0][i]==-1.00){s1=s1-Ad[0][j];}
				else{printf("This is fuck\n");}
			}
			else
			{
				if(M[i][k]==1.00){s1=s1+Ad[k][j];}
				else if(M[i][k]==-1.00){s1=s1-Ad[k][j];}
				else{printf("This is fuck\n");}
			}
			B[i][j]=s1;
		}
		k=k+1;
	}
	//Transpose(n);
	//Transposet(n);
	for(i=0;i<n;i++)
	{
		k=1;
		for(j=0;j<n;j++)
		{
			double s2=0.00;
			if(j!=n-1)
			{
				if(Mt[k][j]==1.00){s2=s2+B[i][k];}
				else if(Mt[k][j]==-1.00){s2=s2-B[i][k];}
				k=k+1;
			}
			else
			{
				if(Mt[0][j]==1.00){s2=s2-B[i][0];}
				else if(Mt[0][j]==-1.00){s2=s2+B[i][0];}
				//k=k+1;
			}
			C[i][j]=s2;
			//putw(C[i][j],fptr);
			printf("%lf\t",C[i][j]);
			fprintf(fptr,"%lf ",C[i][j]);
		}
		printf("\n");
	}
	//Transposet(n);
}
int main()
{
	clock_t start, end;
    	double cpu_time_used;
    	start = clock();
	fptr = fopen("matrix", "w");
	fpt = fopen("keys", "w");
	fp = fopen("ot_transpose","w");
	int z=input();
	if(z!=0)
	{
		keygen();
		encrypt();
	}
	fclose(fptr);
	fclose(fpt);
	fclose(fp);
	end = clock();
      	cpu_time_used = (double) (end - start)/(double)(CLOCKS_PER_SEC);
      	printf("Running time of code is %lf\n\n",cpu_time_used);
	return 0;
}
