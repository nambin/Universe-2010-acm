#include <stdafx.h>
#include <Basis/TestFunction_TBB.h>

#ifdef USE_TBB

namespace basis {

const size_t L = 150;
const size_t M = 225;
const size_t N = 300;

void SerialMatrixMultiply( float c[M][N], float a[M][L], float b[L][N] ) {
	for( size_t i=0; i<M; ++i ) {
		for( size_t j=0; j<N; ++j ) {
			float sum = 0;
			for( size_t k=0; k<L; ++k )
				sum += a[i][k]*b[k][j];
			c[i][j] = sum;
		}
	}
}

class MatrixMultiplyBody2D {
	float (*my_a)[L];
	float (*my_b)[N];
	float (*my_c)[N];
public:
	void operator( )( const blocked_range2d<size_t>& r ) const {
		float (*a)[L] = my_a; // a,b,c used in example to emphasize
		float (*b)[N] = my_b; // commonality with serial code
		float (*c)[N] = my_c;
		for( size_t i=r.rows().begin(); i!=r.rows().end( ); ++i ){
			for( size_t j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				float sum = 0;
				for( size_t k=0; k<L; ++k )
					sum += a[i][k]*b[k][j];
				c[i][j] = sum;
			}
		}
	}
	MatrixMultiplyBody2D( float c[M][N], float a[M][L], float b[L][N] ) : my_a(a), my_b(b), my_c(c)	{}
};

void ParallelMatrixMultiply(float c[M][N], float a[M][L], float b[L][N]){
	int gr1 = 16;
	int gr2 = 32;

	//parallel_for( blocked_range2d<size_t>(0, M, 0, N), MatrixMultiplyBody2D(c,a,b), auto_partitioner() );
	//parallel_for( blocked_range2d<size_t>(0, M, 1, 0, N, 1), MatrixMultiplyBody2D(c,a,b), auto_partitioner() );
	//parallel_for( blocked_range2d<size_t>(0, M, gr1, 0, N, gr2), MatrixMultiplyBody2D(c,a,b), auto_partitioner() );
	//parallel_for( blocked_range2d<size_t>(0, M, gr1, 0, N, gr2), MatrixMultiplyBody2D(c,a,b) );
}

void initializeMatrix(float a[M][L], float b[L][N]) {
	for(int i=0;i<M;++i) {
		for(int j=0;j<L;++j) {
			a[i][j] = ((float)rand() / RAND_MAX);
		}
	}

	for(int i=0;i<L;++i) {
		for(int j=0;j<N;++j) {
			b[i][j] = ((float)rand() / RAND_MAX);
		}
	}	
	cout << " Done building string." << endl;
}

void TestMatrixMultiply() {
	int nThreads = 3;
	int nIter = 150;

	TBBProgramming tbbprogram;
	tbbprogram.initialze(nThreads);
	
	float ret[M][N], retP[M][N], mat0[M][L], mat1[L][N];
	
	tick_count serial_t0 = tick_count::now();
	for(int k=0;k<nIter;++k)
		SerialMatrixMultiply(ret, mat0, mat1);
	tick_count serial_t1 = tick_count::now();
	cout << " Done with serial version." << endl;

	tick_count parallel_t0 = tick_count::now();
	for(int k=0;k<nIter;++k)
		ParallelMatrixMultiply(retP, mat0, mat1);
	tick_count parallel_t1 = tick_count::now();
	cout << " Done with parallel version." << endl;

	for (size_t i = 0; i < M; ++i) {
		for(size_t j=0;j<N;++j) {
			if(ret[i][j] != retP[i][j]) {
				cout << "ERROR: Serial and Parallel Results are Different!" << endl;
			}
		}			
	}
	cout << " Done validating results." << endl;

	cout << "Serial version ran in " << (serial_t1 - serial_t0).seconds() << " seconds" << endl
		 << "Parallel version ran in " <<  (parallel_t1 - parallel_t0).seconds() << " seconds" << endl
		 << "Resulting in a speedup of " << (serial_t1 - serial_t0).seconds() / (parallel_t1 - parallel_t0).seconds() << endl;

	tbbprogram.terminate();
}

};		// basis

#endif		// USE_TBB