#ifndef _MATRIX_VECTOR_OPERATION_H_
#define _MATRIX_VECTOR_OPERATION_H_

// M * a
template <typename MATRIX, typename VECTOR>
struct Vector_M_mul_V {
	MATRIX const & M ;
	VECTOR const & a ;
	Vector_M_mul_V(MATRIX const & _M, VECTOR const & _a) : M(_M), a(_a) {}
} ;

// s * M * a
template <typename SCALAR, typename MATRIX, typename VECTOR>
struct Vector_S_mul_M_mul_V {
	SCALAR const & s ;
	MATRIX const & M ;
	VECTOR const & a ;
	Vector_S_mul_M_mul_V(SCALAR const & _s,	MATRIX const & _M, VECTOR const & _a) : s(_s), M(_M), a(_a) {}
} ;

// a + M * b
template <typename MATRIX, typename VECTOR>
struct Vector_V_sum_M_mul_V {
	VECTOR const & a ;
	MATRIX const & M ;
	VECTOR const & b ;
	Vector_V_sum_M_mul_V(VECTOR const & _a, MATRIX const & _M, VECTOR const & _b) : a(_a), M(_M), b(_b) {}
} ;

// a - M * b
template <typename MATRIX, typename VECTOR>
struct Vector_V_sub_M_mul_V {
	VECTOR const & a ;
	MATRIX const & M ;
	VECTOR const & b ;
	Vector_V_sub_M_mul_V(VECTOR const & _a, MATRIX const & _M, VECTOR const & _b) : a(_a), M(_M), b(_b) {}
} ;

// a + s * M * b
template <typename SCALAR, typename MATRIX, typename VECTOR>
struct Vector_V_sum_S_mul_M_mul_V {
	VECTOR const & a ;
	SCALAR const & s ;
	MATRIX const & M ;
	VECTOR const & b ;
	Vector_V_sum_S_mul_M_mul_V(VECTOR const & _a, SCALAR const & _s, MATRIX const & _M, VECTOR const & _b) : a(_a), s(_s), M(_M), b(_b) {}
} ;

// a - s * M * b
template <typename SCALAR, typename MATRIX, typename VECTOR>
struct Vector_V_sub_S_mul_M_mul_V {
	VECTOR const & a ;
	SCALAR const & s ;
	MATRIX const & M ;
	VECTOR const & b ;
	Vector_V_sub_S_mul_M_mul_V(VECTOR const & _a, SCALAR const & _s, MATRIX const & _M, VECTOR const & _b) : a(_a), s(_s), M(_M), b(_b) {}
} ;

#endif