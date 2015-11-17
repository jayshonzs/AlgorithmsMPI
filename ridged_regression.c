#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int inv(char *srcfile, char *invfile, int n);
int multiply(char *leftfile, char *rightfile, char *targetfile, int left_m, int left_n, int right_m, int right_n);
int transpose(char *srcfile, char *targetfile, int m, int n);
int add(char *file1, char *file2, char *targetfile, int m, int n, int add);

int create_diag(int n, float lambda, char *outfile)
{
	int i = 0, j = 0;
	FILE* fp = fopen(outfile, "wb");
	float val = 0;
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			if(i == j)
			{
				val = lambda;
			}
			else
			{
				val = 0;
			}
			fwrite(&val, sizeof(float), 1, fp);
		}
	}
	
	fclose(fp);
	
	return 0;
}

int ridged_regression(char *xfile, char *yfile, int m, int n, float lambda)
{	
	create_diag(n, lambda, "diag.dat");

	transpose(xfile, "XT.dat", m, n);
	
	multiply("XT.dat", xfile, "XTX.dat", n, m, m, n);

	add("XTX.dat", "diag.dat", "M.dat", n, n, 1);
	
	inv("M.dat", "INV.dat", n);
	
	multiply("INV.dat", "XT.dat", "MXT.dat", n, n, n, m);
	
	multiply("MXT.dat", yfile, "beta.para", n, m, m, 1);
	
	return 0;
}
