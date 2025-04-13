#ifndef _ARITHMETIC_SOLVER_ITERATION_H_
#define _ARITHMETIC_SOLVER_ITERATION_H_

#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../Vector/VectorNorm.h"

namespace arithmetic {

template<typename VEC> 
class Basic_Iteration : public basis::ParallelComputable {
public : 
	typedef typename VEC::value_type value_type;

	Basic_Iteration(size_t iter, value_type eps, const VectorNorm<VEC>* norm = NULL);
	~Basic_Iteration() { if(norm) delete norm; }

	void setComputingOption(const ParallelComputable& rhs);

	const size_t get_num_iter() const { return num_iter; }
	
	const bool finished(const VEC&) const;

	void post_process() const;

public :
	const VectorNorm<VEC>*	norm;

	size_t					max_iter;
	value_type				eps;

	mutable value_type		error;
	mutable size_t			num_iter;
};

template<typename VEC>
Basic_Iteration<VEC>::Basic_Iteration(size_t iter, value_type _eps, const VectorNorm<VEC>* _norm) {
	if(_norm == NULL)	norm = new VectorNorm_Inf<VEC>();
	else				norm = _norm;

	((VectorNorm<VEC>*)norm)->setComputingOption(*this);

	max_iter = iter;
	eps = _eps;

	error = 0;
	num_iter = 0;
}

template<typename VEC>
void Basic_Iteration<VEC>::setComputingOption(const ParallelComputable& rhs) {
	ParallelComputable::setComputingOption(rhs);
	if(norm) ((VectorNorm<VEC>*)norm)->setComputingOption(rhs);
}

template<typename VEC>
const bool Basic_Iteration<VEC>::finished(const VEC& vec) const {
	if(++num_iter == max_iter) 
		return true;

	error = (*norm)(vec);
	return (error <= eps);
}

template<typename VEC> 
void Basic_Iteration<VEC>::post_process() const {
}

};		// arithmetic

#endif