#pragma once

#include <iostream>
#include "../../../Basis/FileIO/FileIO.h"
#include "../../Arithmetic/ArrayUtil/ArrayUtil2.h"
#include "../../Arithmetic/ArrayUtil/ArrayUtil3.h"
#include "../../Arithmetic/MathLibrary/GaussQuadrature.h"
#include "../../Arithmetic/MathLibrary/Interpolation_upwind.h"
#include "../../Arithmetic/GridLibrary/GridLibrary.h"
#include "../MarchingCube/MarchingTriangle.h"
#include "../DistanceFieldLibrary/DistanceFields.h"
#include "../Levelset/Levelset2D.h"
#include "../Levelset/Levelset3D.h"

#define SRL_DEFAULT_QUAD_NODE 5

namespace geometryspace {

template<typename T, int N> class SRL;
template<typename T, int N> class SRL2D;
template<typename T, int N> class SRL3D;
template<typename T, int N> class SRLNode2D;
template<typename T, int N> class SRLNode3D;

template<typename T, int N=SRL_DEFAULT_QUAD_NODE>
class SRL {
public : 
	static const int n_quadnode			= N;
	static const int n_quadnode_1		= N - 1;
	static const int n_quadnode_2		= N - 2;
	static const int n_quadnode_half	= (N - 1) / 2;

	enum { NOT_REFINE = -1, NEED_REFINE = -2 };

	typedef arithmetic::grid_index_type index_type;
	typedef T							value_type;
	
public : 
	SRL() {}

	void release();
	void set(T h, T narrowband, bool uniform_quadrature);

	void get_weight(T* weight, const T x) const;	// x : normalized coordinate
	const int get_sub_index(const T x) const;		// x : normalized coordinate

	virtual const bool exportBinary(std::ofstream&) const;
	virtual const bool importBinary(std::ifstream&);

	virtual const bool exportBinary_Compressed(std::ofstream&) const = 0;
	virtual const bool importBinary_Compressed(std::ifstream&) = 0;

public : 
	bool uniform_quadrature;

	T	narrowband;					// Normalized Value	(0,1,2,...)
	T	narrowbandWidth;			// Unnormalized Value (0,h,2h,...)			

	T	r_sub[N];					// unitized position of sub-cell node, i.e. r_sub[0] = 0; r_sub[n_quadnode-1] = 1
	T	h_norm_sub[N-1];
	T	hInv_norm_sub[N-1];

	T	pos_sub[N];					// pos_sub[i] = r_sub[i] * h 
	T	h_sub[N-1];					// h_sub[i] = pos_sub[i+1] - pos_sub[i]
	T	hInv_sub[N-1];				// hInv_sub[i] = 1 / h_sub[i]	

	T	qInv[N];

	arithmetic::WENO_Adaptive<T>	weno_adaptive[N];
	arithmetic::WENO3rd_Adaptive<T> weno3rd_adaptive[N];

public : 
	class SRL_Interpolator {
	public : 
		typedef typename SRL<T,N>::index_type index_type;

		SRL_Interpolator(bool clip) { clipping = clip; num_itrp = num_clipping = 0; }

		const T operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const;
		const T operator()(const SRL2D<T,N>& srl, index_type i[2], const Vector2<T>& x) const { return (*this)(srl, i[0], i[1], x); }

		const T operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const;
		const T operator()(const SRL3D<T,N>& srl, index_type i[3], const Vector3<T>& x) const { return (*this)(srl, i[0], i[1], i[2], x); }

		bool clipping;

		mutable unsigned int	num_itrp, num_clipping;
	};

	class SRL_Lerp {
	public : 
		typedef typename SRL<T,N>::index_type index_type;

		SRL_Lerp() {}

		const T operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const;
		const T operator()(const SRL2D<T,N>& srl, index_type i[2], const Vector2<T>& x) const { return (*this)(srl, i[0], i[1], x); }

		const T operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const;
		const T operator()(const SRL3D<T,N>& srl, index_type i[3], const Vector3<T>& x) const { return (*this)(srl, i[0], i[1], i[2], x); }

		template<typename T2> const T2 operator()(const SRL3D<T,N>&, const SRLNode3D<T2,N>&, const Vector3<T>& x) const;
	};

	class SRL_Coarse_Lerp {
	public : 
		typedef typename SRL<T,N>::index_type index_type;

		SRL_Coarse_Lerp() {}

		const T operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const;
		const T operator()(const SRL2D<T,N>& srl, index_type i[2], const Vector2<T>& x) const { return (*this)(srl, i[0], i[1], x); }

		const T operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const;
		const T operator()(const SRL3D<T,N>& srl, index_type i[3], const Vector3<T>& x) const { return (*this)(srl, i[0], i[1], i[2], x); }

		template<typename T2> const T2 operator()(const Grid3<T2,T>&, index_type i, index_type j, index_type k, const Vector3<T>& x) const;
	};
};

template<typename T, int N>
void SRL<T,N>::release() {}

template<typename T, int N>
void SRL<T,N>::set(T h, T _narrowband, bool _uniform_quadrature) {
	if(N % 2 != 1 || N != 5) 
		basis::throwError("Error in SRL<T,N>::set(T h, T _narrowband)");

	uniform_quadrature = _uniform_quadrature;

	narrowband = _narrowband;
	narrowbandWidth = narrowband * h;

	// Set-up Quadrature Points	
	r_sub[0] = 0; r_sub[N-1] = 1;
	if(uniform_quadrature == false) {
		double x[N], w[N];		
		gauss_lobatto_legendre(N, x, w, 0.0, 1.0);

		for(int k=1;k<N-1;++k) r_sub[k] = T(x[k]);	
	}
	else {
		for(int k=1;k<N-1;++k) r_sub[k] = T(k) / T(N-1);	
	}	

	for(int i=0;i<N;++i) {
		pos_sub[i] = r_sub[i] * h;
		qInv[i] = 1;

		for(int j=0;j<N;++j) {
			if(i==j) continue;
			qInv[i] /= (r_sub[i] - r_sub[j]);
		}
	}
	for(int i=0;i<N-1;++i) { 
		h_norm_sub[i] = r_sub[i+1] - r_sub[i];	hInv_norm_sub[i] = T(1.0) / h_norm_sub[i]; 
		h_sub[i] = pos_sub[i+1] - pos_sub[i]; 	hInv_sub[i] = T(1.0) / h_sub[i]; 
	}

	T h_array[6];
	for(int i=0;i<N;++i) {
		for(int k=0;k<6;++k) 
			h_array[k] = h_sub[(i - 3 + k + N-1 + N-1) % (N-1)];
		weno_adaptive[i].set(h_array);
	}	

	T h_array2[4];
	for(int i=0;i<N;++i) {
		for(int k=0;k<4;++k) 
			h_array2[k] = h_sub[(i - 2 + k + N-1 + N-1) % (N-1)];
		weno3rd_adaptive[i].set(h_array2);
	}	

	// Auxiliary Test
//	log_info << "SRL R_Sub : "; for(int i=0;i<N;++i)   log_info << r_sub[i] << " "; cout << endl;
//	log_info << "SRL H_Sub : "; for(int i=0;i<N-1;++i) log_info << h_sub[i] << " "; cout << endl;	
}

template<typename T, int N>
void SRL<T,N>::get_weight(T* w, const T x) const {
	const T eps = 1.0e-6;

	T diff[10], diff_mul = 1;
	for(index_type i=0;i<N;++i) { w[i] = 0;											}
	for(index_type i=0;i<N;++i) { if(::abs(x - r_sub[i]) < eps) { w[i] = 1; return; }	}
	for(index_type i=0;i<N;++i) { diff[i] = x - r_sub[i]; diff_mul *= diff[i];			}
	for(index_type i=0;i<N;++i) { w[i] = diff_mul / diff[i] * qInv[i]; }
}

template<typename T, int N>
const int SRL<T,N>::get_sub_index(const T x) const {
	for(int i=1;i<N-1;++i) {
		if(x < r_sub[i]) 
			return i-1;
	}
	return N - 2;
}

template<typename T, int N> 
const bool SRL<T,N>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(&uniform_quadrature), sizeof(bool));

	ofs.write((char*)(&narrowband), sizeof(T));
	ofs.write((char*)(&narrowbandWidth), sizeof(T));

	ofs.write((char*)(r_sub), sizeof(T) * N);
	ofs.write((char*)(h_norm_sub), sizeof(T) * (N-1));
	ofs.write((char*)(hInv_norm_sub), sizeof(T) * (N-1));

	ofs.write((char*)(pos_sub), sizeof(T) * N);
	ofs.write((char*)(h_sub), sizeof(T) * (N-1));
	ofs.write((char*)(hInv_sub), sizeof(T) * (N-1));

	ofs.write((char*)(qInv), sizeof(T) * N);

	for(int i=0;i<N;++i)
		weno_adaptive[i].exportBinaryFile(ofs);

	return true;
}

template<typename T, int N> 
const bool SRL<T,N>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(&uniform_quadrature), sizeof(bool));

	ifs.read((char*)(&narrowband), sizeof(T));
	ifs.read((char*)(&narrowbandWidth), sizeof(T));

	ifs.read((char*)(r_sub), sizeof(T) * N);
	ifs.read((char*)(h_norm_sub), sizeof(T) * (N-1));
	ifs.read((char*)(hInv_norm_sub), sizeof(T) * (N-1));

	ifs.read((char*)(pos_sub), sizeof(T) * N);
	ifs.read((char*)(h_sub), sizeof(T) * (N-1));
	ifs.read((char*)(hInv_sub), sizeof(T) * (N-1));

	ifs.read((char*)(qInv), sizeof(T) * N);

	for(int i=0;i<N;++i)
		weno_adaptive[i].importBinaryFile(ifs);

	return true;
}

};		// geometryspace

#include "SRL_Interpolation.h"
