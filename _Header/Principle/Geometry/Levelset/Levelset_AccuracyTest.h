#pragma once

#include "Levelset2D.h"
#include "Levelset3D.h"
#include "../SRL/SRL2D.h"

namespace geometryspace {

template<typename T> 
class Levelset_AccuracyTest {
public : 
	typedef typename Levelset2D<T>::index_type index_type;

	Levelset_AccuracyTest(T _narrowband = 3) : narrowband(_narrowband) {}
	void set_narrowband(T _narrowband) { narrowband = _narrowband; }

	virtual const T norm(const Levelset2D<T>&, const DistanceField2D<T>&) const { basis::throwError("Not Implemented _ Levelset_AccuracyTest::norm()"); return -1; }
	virtual const T norm(const Levelset3D<T>&, const DistanceField3D<T>&) const { basis::throwError("Not Implemented _ Levelset_AccuracyTest::norm()"); return -1; }

	virtual const T norm(const SRL2D<T>&, const DistanceField2D<T>&) const { basis::throwError("Not Implemented _ Levelset_AccuracyTest::norm()"); return -1; }

	T narrowband;		// normalized value
};

template<typename T> 
class Levelset_InfNorm : public Levelset_AccuracyTest<T> {
public : 
	Levelset_InfNorm(T _narrowband) : Levelset_AccuracyTest(_narrowband) {}

	const T norm(const Levelset2D<T>& lv, const DistanceField2D<T>& field) const {
		T band = narrowband * lv.h;
		T max_abs = 0, diff_abs;
		T exact;

		for(index_type i=1;i<lv.phi.nX-1;++i) {
			for(index_type j=1;j<lv.phi.nY-1;++j) {
				exact = field(lv.phi.pos(i,j));

				if(lv.phi(i,j) > band || lv.phi(i,j) < - band) 
				//if( exact > band || exact < - band) 
					continue;

				diff_abs = ::abs(lv.phi(i,j) - exact);
				if(max_abs < diff_abs) {
					max_abs = diff_abs;
					max_i = i, max_j = j;
					max_error_phi = lv.phi(i,j); max_error_exact = exact;
				}
			}
		}

		return max_abs;
	}

	const T norm(const SRL2D<T>& srl, const DistanceField2D<T>& field) const {
		T band = narrowband * srl.h;
		T max_abs = 0, diff_abs;
		T exact;
		
		for( index_type i = 0; i < srl.phi.nX-1 ; ++i ) {	
			for( index_type j = 0; j < srl.phi.nY-1 ; ++j ) {	
				if(srl.indices(i,j) == SRL<T>::NOT_REFINE) 
					continue;

				const SRLNode2D<T>&	src_node = srl.node(i,j);

				Vector2<unsigned int> rs, re;
				srl.get_updating_subgrid_index(rs, re, src_node.index);

				for(unsigned int p=rs[0];p<re[0];++p) {
					for(unsigned int q=rs[1];q<re[1];++q) {
						exact = field(srl.phi.pos(i,j) + Vector2<real>(srl.pos_sub[p], srl.pos_sub[q]));
						if(src_node(p,q) > band || src_node(p,q) < - band) 
							continue;

						diff_abs = ::abs(src_node(p,q) - exact);
						if(max_abs < diff_abs) {
							max_abs = diff_abs;

							max_i = i, max_j = j, max_p = p, max_q = q;
							max_error_phi = src_node(p,q);
							max_error_exact = exact;
						}
					}
				}
			}
		}
		
		return max_abs;
	}

	mutable index_type		max_i, max_j;
	mutable unsigned int	max_p, max_q;
	mutable T				max_error_phi, max_error_exact;
};

template<typename T> 
class Levelset_L1Norm : public Levelset_AccuracyTest<T> {
public : 
	Levelset_L1Norm(T _narrowband) : Levelset_AccuracyTest(_narrowband) {}

	const T norm(const Levelset2D<T>& lv, const DistanceField2D<T>& field) const {
		T band = narrowband * lv.h;
		T sum_abs = 0, exact;
		unsigned int n_sample = 0;

		for(index_type i=1;i<lv.phi.nX-1;++i) {
			for(index_type j=1;j<lv.phi.nY-1;++j) {
				exact = field(lv.phi.pos(i,j));

				if(lv.phi(i,j) > band || lv.phi(i,j) < - band) 
				//if(exact > band || exact < - band) 
					continue;
				
				sum_abs += ::abs(lv.phi(i,j) - exact);
				n_sample++;
			}
		}

		return sum_abs / n_sample;
	}

	const T norm(const SRL2D<T>& srl, const DistanceField2D<T>& field) const {
		T band = narrowband * srl.h;
		T sum_abs = 0, exact;
		unsigned int n_sample = 0;

		for( index_type i = 0; i < srl.phi.nX-1 ; ++i ) {	
			for( index_type j = 0; j < srl.phi.nY-1 ; ++j ) {	
				if(srl.indices(i,j) == SRL<T>::NOT_REFINE) 
					continue;

				const SRLNode2D<T>&	src_node = srl.node(i,j);

				Vector2<unsigned int> rs, re;
				srl.get_updating_subgrid_index(rs, re, src_node.index);

				for(unsigned int p=rs[0];p<re[0];++p) {
					for(unsigned int q=rs[1];q<re[1];++q) {
						exact = field(srl.phi.pos(i,j) + Vector2<real>(srl.pos_sub[p], srl.pos_sub[q]));
						if(src_node(p,q) > band || src_node(p,q) < - band) 
							continue;

						sum_abs += ::abs(src_node(p,q) - exact);
						n_sample++;
					}
				}
			}
		}


		return sum_abs / n_sample;
	}
};

};		// geometryspace