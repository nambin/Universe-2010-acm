#ifndef _MATRIX_OPERATION_H_
#define _MATRIX_OPERATION_H_

namespace arithmetic {

template<typename MATRIX, typename VECX, typename VECY>
void mul(VECY& y, const MATRIX& A, const VECX& x) {
	A.mul(y, x);	// y = Ax
}

// -M
template <typename MATRIX>
struct Matrix_neg_M {
	MATRIX const & M ;
	Matrix_neg_M(MATRIX const & _M) : M(_M) {}
} ;

// M^t
template <typename MATRIX>
struct Matrix_trans_M {
	MATRIX const & M ;
	Matrix_trans_M(MATRIX const & _M) : M(_M) {}
} ;

// M + N
template <typename MATRIX>
struct Matrix_M_add_M {
	MATRIX const & M1 ;
	MATRIX const & M2 ;
	Matrix_M_add_M(MATRIX const & _M1, MATRIX const & _M2) : M1(_M1), M2(_M2) {}
} ;

// M - N
template <typename MATRIX>
struct Matrix_M_dif_M {
	MATRIX const & M1 ;
	MATRIX const & M2 ;
	Matrix_M_dif_M(MATRIX const & _M1, MATRIX const & _M2) : M1(_M1), M2(_M2) {}
} ;

// M * s
template <typename MATRIX, typename SCALAR>
struct Matrix_M_mul_S {
	MATRIX const & M ;
	SCALAR const & s ;
	Matrix_M_mul_S(MATRIX const & _M, SCALAR const & _s) : M(_M), s(_s) {}
} ;

// M * N
template <typename MATRIX>
struct Matrix_M_mul_M {
	MATRIX const & M1 ;
	MATRIX const & M2 ;
	Matrix_M_mul_M(MATRIX const & _M1, MATRIX const & _M2) : M1(_M1), M2(_M2) {}
};

};		// arithmetic

#endif