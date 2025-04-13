#pragma once

#include "../MacGrid/MacGrid2D.h"
#include "../MacGrid/MacGridLibrary.h"
#include "SRL.h"

namespace geometryspace {

template<typename T, int N=SRL_DEFAULT_QUAD_NODE>
class SRLNode2D {
	typedef typename SRL<T,N>::index_type index_type;	

public : 
	SRLNode2D() { index[0] = index[1] = -1; }
	SRLNode2D(index_type ii, index_type jj) { index[0] = ii; index[1] = jj; }

	inline const T	operator()(int i, int j) const		{ return val[i][j]; }
	inline T&		operator()(int i, int j)			{ return val[i][j]; }

	inline const T	operator()(const int i[2]) const	{ return val[i[0]][i[1]]; }
	inline T&		operator()(const int i[2])			{ return val[i[0]][i[1]]; }

	inline const T	operator()(unsigned int i, unsigned int j) const	{ return val[i][j]; }
	inline T&		operator()(unsigned int i, unsigned int j)			{ return val[i][j]; }

	inline const T	operator()(const unsigned int i[2]) const			{ return val[i[0]][i[1]]; }
	inline T&		operator()(const unsigned int i[2])					{ return val[i[0]][i[1]]; }

	void set(index_type ii, index_type jj)	{ index[0] = ii; index[1] = jj;								}
	void set(const T v)						{ for(int i=0;i<N;++i) for(int j=0;j<N;++j) val[i][j] = v;	}

	void neumann_x_lower();
	void neumann_x_upper();
	void neumann_y_lower();
	void neumann_y_upper();

	const bool has_coarse_interface() const;
	const bool has_interface() const;	

	const bool exportBinary(std::ofstream&) const;
	const bool importBinary(std::ifstream&);

public : 
	index_type			index[2];
	T 					val[N][N];	
};

template<typename T, int N>
void SRLNode2D<T,N>::neumann_x_lower() {
	int half = N/2;
	for(int i=0;i<half;++i) 
		for(int j=0;j<N;++j) 
			val[i][j] = val[half][j];
}

template<typename T, int N>
void SRLNode2D<T,N>::neumann_x_upper() {
	int half = N/2;
	for(int i=half+1;i<N;++i) 
		for(int j=0;j<N;++j) 
			val[i][j] = val[half][j];
}

template<typename T, int N>
void SRLNode2D<T,N>::neumann_y_lower() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=0;j<half;++j) 
			val[i][j] = val[i][half];
}

template<typename T, int N>
void SRLNode2D<T,N>::neumann_y_upper() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=half+1;j<N;++j) 
			val[i][j] = val[i][half];

}

template<typename T, int N>
const bool SRLNode2D<T,N>::has_coarse_interface() const {
	if( val[0][0] >  0 && val[0][N-1] >  0 && val[N-1][0] >  0 && val[N-1][N-1] >  0 )	return false;
	if( val[0][0] <= 0 && val[0][N-1] <= 0 && val[N-1][0] <= 0 && val[N-1][N-1] <= 0 )	return false;

	return true;
}

template<typename T, int N>
const bool SRLNode2D<T,N>::has_interface() const {
	if(val[0][0] > 0) {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				if(val[i][j] <= 0)
					return true;
			}
		}
		return false;
	}
	else {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				if(val[i][j] > 0)
					return true;
			}
		}
		return false;
	}
}

template<typename T, int N=SRL_DEFAULT_QUAD_NODE>
class SRLNode3D {
	typedef typename SRL<T,N>::index_type index_type;	

public : 
	SRLNode3D() { index[0] = index[1] = index[2] = -1; }
	SRLNode3D(index_type ii, index_type jj, index_type kk) { index[0] = ii; index[1] = jj; index[2] = kk; }

	inline const T	operator()(int i, int j, int k) const	{ return val[i][j][k]; }
	inline T&		operator()(int i, int j, int k)			{ return val[i][j][k]; }

	inline const T	operator()(const int i[3]) const	{ return val[i[0]][i[1]][i[2]]; }
	inline T&		operator()(const int i[3])			{ return val[i[0]][i[1]][i[2]]; }

	inline const T	operator()(unsigned int i, unsigned int j, unsigned int k) const	{ return val[i][j][k]; }
	inline T&		operator()(unsigned int i, unsigned int j, unsigned int k)			{ return val[i][j][k]; }

	inline const T	operator()(const unsigned int i[3]) const	{ return val[i[0]][i[1]][i[2]]; }
	inline T&		operator()(const unsigned int i[3])			{ return val[i[0]][i[1]][i[2]]; }

	void set(index_type ii, index_type jj, index_type kk)	{ index[0] = ii; index[1] = jj; index[2] = kk; }
	void set(const T v)										{ for(int i=0;i<N;++i) for(int j=0;j<N;++j) for(int k=0;k<N;++k) val[i][j][k] = v;	}

	void neumann_x_lower();
	void neumann_x_upper();
	void neumann_y_lower();
	void neumann_y_upper();
	void neumann_z_lower();
	void neumann_z_upper();

	const bool has_coarse_interface() const;
	const bool has_interface() const;	

	const bool exportBinary(std::ofstream&) const;
	const bool importBinary(std::ifstream&);

public : 
	index_type			index[3];
	T 					val[N][N][N];	
};

template<typename T, int N>
void SRLNode3D<T,N>::neumann_x_lower() {
	int half = N/2;
	for(int i=0;i<half;++i) 
		for(int j=0;j<N;++j) 
			for(int k=0;k<N;++k)
				val[i][j][k] = val[half][j][k];
}

template<typename T, int N>
void SRLNode3D<T,N>::neumann_x_upper() {
	int half = N/2;
	for(int i=half+1;i<N;++i) 
		for(int j=0;j<N;++j) 
			for(int k=0;k<N;++k) 
				val[i][j][k] = val[half][j][k];
}

template<typename T, int N>
void SRLNode3D<T,N>::neumann_y_lower() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=0;j<half;++j) 
			for(int k=0;k<N;++k) 
				val[i][j][k] = val[i][half][k];
}

template<typename T, int N>
void SRLNode3D<T,N>::neumann_y_upper() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=half+1;j<N;++j) 
			for(int k=0;k<N;++k) 
				val[i][j][k] = val[i][half][k];
}

template<typename T, int N>
void SRLNode3D<T,N>::neumann_z_lower() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=0;j<N;++j) 
			for(int k=0;k<half;++k) 
				val[i][j][k] = val[i][j][half];
}

template<typename T, int N>
void SRLNode3D<T,N>::neumann_z_upper() {
	int half = N/2;
	for(int i=0;i<N;++i) 
		for(int j=0;j<N;++j) 
			for(int k=half+1;k<N;++k) 
				val[i][j][k] = val[i][j][half];
}

template<typename T, int N>
const bool SRLNode3D<T,N>::has_coarse_interface() const {
	if( val[0][0][0]   >  0 && val[0][N-1][0]   >  0 && val[N-1][0][0]   >  0 && val[N-1][N-1][0]   >  0 
	 && val[0][0][N-1] >  0 && val[0][N-1][N-1] >  0 && val[N-1][0][N-1] >  0 && val[N-1][N-1][N-1] >  0 )		return false;

	if( val[0][0][0]   <=  0 && val[0][N-1][0]   <=  0 && val[N-1][0][0]   <=  0 && val[N-1][N-1][0]   <=  0 
	 && val[0][0][N-1] <=  0 && val[0][N-1][N-1] <=  0 && val[N-1][0][N-1] <=  0 && val[N-1][N-1][N-1] <=  0 )	return false;

	return true;
}

template<typename T, int N>
const bool SRLNode3D<T,N>::has_interface() const {
	if(val[0][0][0] > 0) {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				for(int k=0;k<N;++k) {
					if(val[i][j][k] <= 0)
						return true;
				}
			}
		}
		return false;
	}
	else {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				for(int k=0;k<N;++k) {
					if(val[i][j][k] > 0)
						return true;
				}
			}
		}
		return false;
	}
}

template<typename T, int N> 
const bool SRLNode2D<T,N>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(&index), sizeof(index_type) * 2);
	ofs.write((char*)(&val), sizeof(T) * N * N);	
	return true;
}

template<typename T, int N> 
const bool SRLNode2D<T,N>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(&index), sizeof(index_type) * 2);
	ifs.read((char*)(&val), sizeof(T) * N * N);	
	return true;
}

template<typename T, int N> 
const bool SRLNode3D<T,N>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(&index), sizeof(index_type) * 3);
	ofs.write((char*)(&val), sizeof(T) * N * N * N);	
	return true;
}

template<typename T, int N> 
const bool SRLNode3D<T,N>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(&index), sizeof(index_type) * 3);
	ifs.read((char*)(&val), sizeof(T) * N * N * N);	
	return true;
}

};		// namespace geometryspace

