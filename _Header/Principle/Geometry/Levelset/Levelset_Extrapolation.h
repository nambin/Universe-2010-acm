#ifndef _GEOMETRYSPACE_LEVELSET_EXTRAPOLATION_H_
#define _GEOMETRYSPACE_LEVELSET_EXTRAPOLATION_H_

#include "Levelset2D.h"
#include "Levelset3D.h"
#include "FMMHeap.h"

namespace geometryspace {

template<typename T>
class Levelset_Extrapolation {
public : 
	Levelset_Extrapolation(T band) : narrowBandWidth(band) {}

	virtual void post_process() const { log_info << "Levelset_Extrapolation _ Have no Post Process" << std::endl; }

	virtual void extrapolate(Grid2<T,T>&, typename MacGrid2D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid2<T,T>&) const = 0;
	virtual void extrapolate(Grid3<T,T>&, typename MacGrid3D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid3<T,T>&) const = 0;	

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ Narrowband : " << narrowBandWidth << endl; }
	static const std::string XmlElementName() { return std::string("Levelset_Extrapolation"); }

	T	narrowBandWidth;		// Normalized Value	(0,1,2,...)
};

template<typename T>
class Levelset_ExtrapolationFMM : public Levelset_Extrapolation<T> {
public : 
	Levelset_ExtrapolationFMM(T band) : Levelset_Extrapolation(band), heap(NULL) {}
	~Levelset_ExtrapolationFMM() { if(heap) delete heap; }

	void extrapolate(Grid2<T,T>&, typename MacGrid2D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid2<T,T>&) const;
	void extrapolate(Grid3<T,T>&, typename MacGrid3D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid3<T,T>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Extrapolation::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_ExtrapolationFMM"); }

private : 
	void extrapolate(Grid2<T,T>&, const GridBoundaryCondition<T,T>&, const Grid2<T,T>& phi) const;
	void extrapolate(Grid3<T,T>&, const GridBoundaryCondition<T,T>&, const Grid3<T,T>& phi) const;
	
	//////////////////////////////////////////////////////////////////////////
	// Variable for Fast Marching
	mutable FMMHeap<T>*			heap;
	mutable grid_index_type		minHeapSize;
};

template<typename T>
class Levelset_Extrapolation_PDE : public Levelset_Extrapolation<T> {
public : 
	Levelset_Extrapolation_PDE(T band, T _cfl = 0.5f) : Levelset_Extrapolation(band), cfl(_cfl) { assert(cfl <= T(1.0)); cfl = MIN(cfl, T(0.5)); }

	void extrapolate(Grid2<T,T>&, typename MacGrid2D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid2<T,T>&) const;
	void extrapolate(Grid3<T,T>&, typename MacGrid3D<T>::FIELD_TYPE type, const GridBoundaryCondition<T,T>&, const Grid3<T,T>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ CFL : " << cfl << endl; Levelset_Extrapolation::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_ExtrapolationPDE"); }

private : 
	void extrapolate(Grid2<T,T>&, const GridBoundaryCondition<T,T>&, const Grid2<T,T>& phi) const;
	void extrapolate(Grid3<T,T>&, const GridBoundaryCondition<T,T>&, const Grid3<T,T>& phi) const;

	T cfl;
};

};		// geometryspace

#include "Levelset_ExtrapolationFMM2D.h"
#include "Levelset_ExtrapolationFMM3D.h"
#include "Levelset_Extrapolation_PDE.h"

#endif