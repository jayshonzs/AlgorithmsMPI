#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

int km_loaddata(double data[300][2], const char* filename)
{
	int ret = 0;
	int i = 0;
	FILE* fp = fopen(filename, "r");
	char buf[256];
	char *res = NULL;
	
	while(1)
	{
		res = fgets(buf, 256, fp);
		if(res == NULL)
		{
			//printf("km_loaddata fgets error!\n");
			break;
		}
		sscanf(buf, "%lf %lf", &data[i][0], &data[i][1]);
		i++;
	}
	
	printf("load %lf %lf\n", data[0][0], data[0][1]);
	printf("load %lf %lf\n", data[299][0], data[299][1]);
	
	fclose(fp);
	
	return ret;
}

int gprint = 0;

double euclidean(double *x1, double *x2, int D)
{
	double ret = 0;
	int d = 0;
	double sum = 0.0;
	
	for(d = 0; d < D; d++)
	{
		sum += pow(x1[d]-x2[d], 2);
	}
	
	ret = sqrt(sum);
	
	if(gprint)
	{
		printf("sum:%lf\n", sum);
		printf("ret:%lf\n", ret);
	}
	
	return ret;
}

int scatter_orig_data(double sbuf[300][2], double rb[75][2], int rank, int numtasks, MPI_Datatype *coordtype)
{
	int src = 0;
	int scnt = 300/numtasks;
	int rcnt = 300/numtasks;
	int error_code;
	char error_string[1024];
	int resultlen;
	
	rb[0][0] = 0;

	//printf("Proc %d, address %llx\n", rank, rb);
	
	error_code = MPI_Scatter(sbuf,		//散发数据缓冲区
							scnt,		//散发数据个数
							*coordtype,	//散发数据类型
							rb,			//接收数据缓冲区
							rcnt,		//接收数据个数
							*coordtype,	//接收数据类型
							src,		//源进程号
							MPI_COMM_WORLD);
	
	if(error_code != MPI_SUCCESS)
	{
		MPI_Error_string(error_code, error_string, &resultlen);
		printf("Proc %d error: %s\n", rank, error_string);
		exit(-1);
	}
	
	printf("Proc %d rb[0][0]=%lf, rb[0][1]=%lf\n", rank, rb[0][0], rb[0][1]);
	printf("Proc %d rb[%d][0]=%lf, rb[%d][1]=%lf\n", rank, rcnt-1, rb[rcnt-1][0], rcnt-1, rb[rcnt-1][1]);
	
	MPI_Barrier(MPI_COMM_WORLD);

	return 0;
}

#if 1
void mysum1(double c_in[3][2], double c_inout[3][2], int *len, MPI_Datatype *dptr)
{
	int i;
	//double c[2];
	
	for(i = 0; i < *len; i++)
	{
		c_inout[i][0] = c_in[i][0] + c_inout[i][0];
		c_inout[i][1] = c_in[i][1] + c_inout[i][1];
	}
}

void mysum2(double c_in[3], double c_inout[3], int *len, MPI_Datatype *dptr)
{
	int i;
	
	for(i = 0; i < *len; i++)
	{
		c_inout[i] = c_in[i] + c_inout[i];
	}
}
#endif

int run_kmeans(double orig_data[300][2], double local_data[75][2], int N, int K, int rank, MPI_Datatype coordtype, double out_centers[3][2])
{
	double centers[3][2], new_centers[3][2], recv_centers[3][2];
	int cluster_sizes[3], recv_sizes[3];
	int i = 0, k = 0;
	double dist = 0, min_dist = 0;
	int best_k = 0;
	int belong[N];
	int temp = 0;
	int indexes[3];
	MPI_Op myOp1, myOp2;
	double differ = 0, res;
	int loop = 0;
	
	MPI_Op_create((MPI_User_function*)mysum1, 1, &myOp1);
	MPI_Op_create((MPI_User_function*)mysum2, 1, &myOp2);
	
	//initialization of centers
	if(rank == 0)
	{
		srand((int)time(0));
		for(k = 0; k < K; k++)
		{
			while(1)
			{
				temp = rand()%N;
				for(i = 0; i < k; i++)
				{
					if(temp == indexes[i])
					{
						break;
					}
				}
				if(i == k)
				{
					break;
				}
			}
			centers[k][0] = orig_data[temp][0];
			centers[k][1] = orig_data[temp][1];
			
			printf("Proc %d initial center[%d][0]=%lf, center[%d][1]=%f\n", rank, k, centers[k][0], k, centers[k][1]);
		}
	}
	
	while(1)
	{
		loop++;
		//printf("loop:%d, rank:%d 0\n",loop, rank);
		MPI_Bcast(centers, K, coordtype, 0/*root process*/, MPI_COMM_WORLD);
		
		if(rank == 1)
		{
			for(k = 0; k < K; k++)
			{
				printf("Proc %d Bcast center[%d][0]=%lf, center[%d][1]=%f\n", rank, k, centers[k][0], k, centers[k][1]);
			}
		}
				
		for(k = 0; k < K; k++)
		{
			cluster_sizes[k] = 0;
			new_centers[k][0] = 0;
			new_centers[k][1] = 0;
		}
		
		for(i = 0; i < N; i++)
		{
			min_dist = 999999999;
			best_k = -1;
			for(k = 0; k < K; k++)
			{
				dist = euclidean(local_data[i], centers[k], 2);
				
				if(dist < min_dist)
				{
					min_dist = dist;
					best_k = k;
				}
			}
			belong[i] = best_k;
		}
		
		for(i = 0; i < N; i++)
		{
			k = belong[i];
			new_centers[k][0] += local_data[i][0];
			new_centers[k][1] += local_data[i][1];
			cluster_sizes[k]++;
		}
		
		/*
		for(k = 0; k < K; k++)
		{
			if(cluster_sizes[k] == 0)
			{
				continue;
			}
			new_centers[k][0] /= cluster_sizes[k];
			new_centers[k][1] /= cluster_sizes[k];
		}
		*/
		
		#if 0
		if(rank == 0)
		{
			for(k = 0; k < K; k++)
			{
				printf("Proc %d before Reduce new_center[%d][0]=%lf, new_center[%d][1]=%f\n", rank, k, new_centers[k][0], k, new_centers[k][1]);
			}
		}
		#endif
		
		MPI_Reduce(new_centers, recv_centers, K, coordtype, myOp1, 0, MPI_COMM_WORLD);
		MPI_Reduce(cluster_sizes, recv_sizes, K, MPI_DOUBLE, myOp2, 0, MPI_COMM_WORLD);
		
		#if 1
		if(rank == 0)
		{
			for(k = 0; k < K; k++)
			{
				printf("Proc %d after Reduce recv_center[%d][0]=%lf, recv_center[%d][1]=%f\n", rank, k, recv_centers[k][0], k, recv_centers[k][1]);
			}
		}
		#endif
		
		if(rank == 0)
		{
			for(k = 0; k < K; k++)
			{
				recv_centers[k][0] /= recv_sizes[k];
				recv_centers[k][1] /= recv_sizes[k];
			}
			
			differ = 0;
			gprint = 1;
			for(k = 0; k < K; k++)
			{
				res = euclidean(centers[k], recv_centers[k], 2);
				differ += res;
				printf("loop: %d, centers[%d][0]=%lf, centers[%d][1]=%lf, recv_centers[%d][0]=%lf, recv_centers[%d][1]=%lf, res=%lf, differ=%lf\n",
						loop, k, centers[k][0], k, centers[k][1], k, recv_centers[k][0], k, recv_centers[k][1], res, differ);
			}
			gprint = 0;
			
			for(k = 0; k < K; k++)
			{
				centers[k][0] = recv_centers[k][0];
				centers[k][1] = recv_centers[k][1];
			}
		}
		
		MPI_Bcast(&differ, 1, MPI_DOUBLE, 0/*root process*/, MPI_COMM_WORLD);
		if(differ < 0.000001)
		{
			break;
		}
	}
	
	for(k = 0; k < K; k++)
	{
		out_centers[k][0] = centers[k][0];
		out_centers[k][1] = centers[k][1];
	}
	
	MPI_Op_free(&myOp1);
	MPI_Op_free(&myOp2);
	
	return 0;
}

int save_result(double orig_data[300][2], double centers[3][2], const char* filename)
{
	int i, k;
	FILE* fp = fopen(filename, "w");
	int best_k = -1;
	double min_dist = 0, dist = 0;
	
	//printf("############ %lf %lf %d\n", orig_data[0][0], orig_data[0][1], best_k);
	
	for(i = 0; i < 300; i++)
	{
		best_k = -1;
		min_dist = 999999999;
		for(k = 0; k < 3; k++)
		{
			dist = euclidean(orig_data[i], centers[k], 2);
			printf("%lf, %lf, %lf, %lf, %lf\n", orig_data[i][0], orig_data[i][1], centers[k][0], centers[k][1], dist);
			if(dist < min_dist)
			{
				min_dist = dist;
				best_k = k;
			}
		}
		fprintf(fp, "%lf %lf %d\n", orig_data[i][0], orig_data[i][1], best_k);
		//printf("%lf %lf %d\n", orig_data[i][0], orig_data[i][1], best_k);
	}
	
	fclose(fp);
	
	return 0;
}

int kmeans(int K, double centers[3][2])
{
	MPI_Datatype coordtype;
	int rank, numtasks;
	double orig_data[300][2];
	double recv_buf[75][2];
	int points_per_proc = 0;
	int i = 0;
	
	//printf("1\n");
	
	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	
	points_per_proc = 300/numtasks;
	
	//printf("2\n");
	
	/*
	orig_data = (double**)malloc(300*sizeof(double*));
	for(i = 0; i < 300; i++)
	{
		orig_data[i] = (double*)malloc(2*sizeof(double));
	}
	
	recv_buf = (double**)malloc(points_per_proc*sizeof(double*));
	for(i = 0; i < points_per_proc; i++)
	{
		recv_buf[i] = (double*)malloc(2*sizeof(double));
	}
	*/
	
	//printf("3\n");
	
	MPI_Type_contiguous(2, MPI_DOUBLE, &coordtype);
	MPI_Type_commit(&coordtype);
	
	//printf("4\n");
	
	//从文件加载数据
	if(rank == 0)
	{
		km_loaddata(orig_data, "gaussian_mix.dat");
	}
	
	//发散数据
	scatter_orig_data(orig_data, recv_buf, rank, numtasks, &coordtype);
	
	//do job
	run_kmeans(orig_data, recv_buf, points_per_proc, K, rank, coordtype, centers);
	
	//done
	if(rank == 0)
	{
		save_result(orig_data, centers, "result.dat");
	}
	
	MPI_Finalize();
	
	return 0;
}
