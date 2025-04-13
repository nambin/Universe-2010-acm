#ifndef _GEOMETRYSPACE_LEVELSET_3D_H_
#define _GEOMETRYSPACE_LEVELSET_3D_H_

#include <vector>
#include "../MacGrid/MacGrid3D.h"
#include "../MarchingCube/MarchingCube.h"
#include "../Surface/TriMesh_Simple.h"
#include "../DistanceFieldLibrary/DistanceFields.h"
#include "../../../Principle/VirtualSpace/Material/Material.h"

namespace geometryspace {

using std::vector;
using geometryspace::TriMesh_Simple;

template<typename T>
class Levelset3D : public MacGrid3D<T> {
public : 
	typedef typename geometryspace::MacGrid3D<T>::index_type index_type;

public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Levelset3D();

	//////////////////////////////////////////////////////////////////////////	
	// Initialization
	virtual void set(index_type nX, index_type nY, index_type nZ, T h);	
	virtual void setPhi(const Grid3<T,T>&);
	virtual void setPhi(const DistanceField3D<T>&);

	virtual void minPhi(const DistanceField3D<T>&);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void createRenderingPrimitive();	
	void createRenderingPrimitive(const arithmetic::Grid3<bool,T>& isbdry);

	const T getTotalVolume() const;
	const T getTotalVolume(const Grid3<bool,T>& isbdry) const;
	const T getVolumeFraction(index_type i, index_type j, index_type k) const { return MacGridLibrary<T>::heaviside(phi(i,j,k) * hInv, 0.0f, 1.0f); }

	void getPhi(Grid3<T,T>& _phi, typename MacGrid3D<T>::FIELD_TYPE type) const { getPhi(_phi, phi, type); }
	static void getPhi(Grid3<T,T>&, const Grid3<T,T>&, typename MacGrid3D<T>::FIELD_TYPE type);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawSurfaceMesh(const Material& mat) const;	

	static void marchingCube(TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>& mesh, const Grid3<T,T>&, const Grid3<bool,T>&);

public : 
	Grid3<T,T>	phi;

	//////////////////////////////////////////////////////////////////////////
	// Variable for Visualization	
	mutable TriMesh_Simple<Vector3<T>, vector<Vector3<T>>> mesh;	
};

template<typename T>
Levelset3D<T>::Levelset3D() {}

template<typename T>
void Levelset3D<T>::set(index_type nX, index_type nY, index_type nZ, T h) {
	MacGrid3D::set(nX, nY, nZ, h);

	const Vector3<T> sOffset(T(-0.5), T(-0.5), T(-0.5));

	phi.init(nX+2, nY+2, nZ+2, h, h, h, sOffset);
	mesh.init();	
}

template<typename T> 
void Levelset3D<T>::setPhi(const Grid3<T,T>& _phi) {
	phi = _phi;
}

template<typename T>
void Levelset3D<T>::setPhi(const DistanceField3D<T>& field) {
	Timer timer;

	int numUpdates = 0;

#ifdef USE_OPENMP
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = phi.nX-1;
	index_type i;
#endif
	for(i = 1; i < my_end; ++i) {	
		Vector3<T> pos;
		for(index_type j=1;j<phi.nY-1;++j)	{
			for(index_type k=1;k<phi.nZ-1;++k)	{
				pos = phi.pos(i,j,k);
				phi(i,j,k) = field.getPhi(pos);
			}			
		}
	}

	phi.setNeumannBoundary();

	timer.logTime("Levelset3D<T>::setPhi()", log_info);
	log_info << "Number of Updates : " << numUpdates << std::endl;
}

template<typename T>
void Levelset3D<T>::minPhi(const DistanceField3D<T>& field) {

#ifdef USE_OPENMP
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = phi.nX-1;
	index_type i;
#endif
	for(i = 1; i < my_end; ++i) {	
		Vector3<T> pos;
		for(index_type j=1;j<phi.nY-1;++j)	{
			for(index_type k=1;k<phi.nZ-1;++k)	{
				pos = phi.pos(i,j,k);
				
				const T _phi = field.getPhi(pos);
				phi(i,j,k) = (phi(i,j,k) < _phi) ? phi(i,j,k) : _phi;
			}			
		}
	}

	phi.setNeumannBoundary();
}

};		// namespace geometryspace

#include "Levelset3D_Utility.h"
#include "Levelset3D_glDraw.h"

#endif

