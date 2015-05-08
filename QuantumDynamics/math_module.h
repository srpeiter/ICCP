#ifndef MATH_MODULE_H
#define MATH_MODULE_H

#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<superlu/slu_zdefs.h>
#include<superlu/slu_dcomplex.h>
#include<cmath>
#include<string.h>

extern "C"
{
	int zgemv_(char* TRANS,int* M,int* N,doublecomplex* ALPHA,doublecomplex* A,int* LDA,doublecomplex* X, int*INCX,
          doublecomplex*  BETA, doublecomplex* Y, int* INCY );

	int zgtsv_(int* N, int* NRHS, doublecomplex* DL, doublecomplex* D, doublecomplex* DU, doublecomplex* B, int* LDB, int* info);

}

namespace math_module
{

template<typename Type>
inline void maketoeplitz(Type** mat, const Type* a, const Type* b, const int N)
{
// making toeplitz matrix
for (int j = 0; j < N; j++)
	for (int k = j; k < N; k++)
		{

		mat[j][k]=a[k-j];
		mat[k][j]=b[k-j];
}

}


template<typename Type>
inline void mat_to_vec( Type* Fortvec, Type** mat,const int N)
{
for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	Fortvec[i*N + j]= mat[j][i];	// transform to column wise for Fortran



}

template<typename Type>
inline void CCS_scheme( Type* data, int* row_ind, int* col_ptr, Type** mat,const int N )
{
	// we know we are dealing with a triangular matrix, but we are implementing a general algortihm,
	// in case we want to implement a 2D case.
	Type* fortvec = new Type[N*N];
	mat_to_vec<Type>(fortvec, mat, N);

	int nnz = 3*N - 2 ; 	// for triangular banded matrices ; number of nonzero elements;
	int j=0;
	int k=0;
	for(int i = 0; i < (N*N) ; i++)
{

	if(fortvec[i].r != 0 | fortvec[i].i !=0)
		{data[j++]=fortvec[i];
		row_ind[k++]= i % N;
		}
}

// we know we are dealing with triangular matrix;
col_ptr[0]=0;
col_ptr[1]=2;

for (int i = 2 ; i < N; i++)
	col_ptr[i] = col_ptr[i-1] + 3;

col_ptr[N]= nnz;

delete [] fortvec;

}


template<typename Type>
inline void MakeComplexSupermatrix(SuperMatrix& A, Type* data, int*row_ind, int* col_ptr, Type** mat, const int N)
{
// setting dimensions
int nnz = 3*N -2;
int nc, nr;
nc =nr = N;

// setting Compressed column storage
CCS_scheme<Type>(  data, row_ind, col_ptr,  mat, N );

zCreate_CompCol_Matrix(&A, nc, nr, nnz, data, row_ind, col_ptr, SLU_NC, SLU_Z, SLU_GE);

}

template<typename Type>
inline void Tofortran(Type* vec, Type** mat, const int N)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{vec[i*N + j].r= mat[j][i].r;	// transform to column wise for Fortran
		vec[i*N + j].i= mat[j][i].i ;}


	}


template<typename Type>
inline void matvec_multiply(Type* out, Type** mat, Type* vec  ,int N)
{
	memset(out , 0, N*sizeof(out[0]));
	Type* mat_temp  = new Type[N*N] ();
	Tofortran<Type> (mat_temp, mat, N);
	char trans = 'N';
	int M = N;
	doublecomplex alpha ; alpha.r=1.0; alpha.i=0;
	int LDA = N;
	int incx = 1;
	doublecomplex beta;  beta.r = 0 ; beta.i = 0;

	int incy = 1;

	zgemv_(&trans,&M,&N,&alpha,mat_temp,&LDA,vec, &incx ,&beta,out,&incy);

	delete [] mat_temp;
}

template<typename Type>
inline void trans_pose(Type** out, Type** mat, const int N)
{
	for(int i = 0 ; i < N; i++)
		for(int j=0; j < N; j++)
			out[i][j]=mat[j][i];
}


template <typename Type>
inline void print(Type** mat, const int N)
{
for (int i = 0; i < N; ++i)
{
  for (int j = 0; j < N; ++j)
	{
		std::cout << std::setw(5) << std::setprecision(4)<<  mat[i][j].r << " + i"<< mat[i][j].i <<"  ";


	}
	std::cout << std::endl;


	    }
std::cout << std::endl;

}

template<typename Type>
inline void print(Type* Fortvec, const int N)
{
	for (int i=0; i < N; i++)
		std::cout << std::setw(3) << Fortvec[i].r  << " + i" << Fortvec[i].i << ' ' ;

	std::cout << std::endl;
}

}





#endif 
