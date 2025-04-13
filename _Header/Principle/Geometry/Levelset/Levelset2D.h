#ifndef _GEOMETRYSPACE_LEVELSET_2D_H_
#define _GEOMETRYSPACE_LEVELSET_2D_H_

#include "../../Arithmetic/Vector/Vector2.h"
#include "../MacGrid/MacGrid2D.h"
#include "../MacGrid/MacGridLibrary.h"
#include "../MarchingCube/MarchingTriangle.h"
#include "../DistanceFieldLibrary/DistanceFields.h"
#include "FMMheap.h"

namespace geometryspace {

using std::vector;

template<typename T>
class Levelset2D : public MacGrid2D<T> {
public : 
	typedef typename geometryspace::MacGrid2D<T>::index_type index_type;

public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Levelset2D();

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	virtual void set(index_type nX, index_type nY, T h);
	virtual void setPhi(const Grid2<T,T>& phi);
	virtual void setPhi(const DistanceField2D<T>&);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void createRenderingPrimitive();
	void createRenderingPrimitive(const arithmetic::Grid2<bool,T>& isbdry);
	
	const T getTotalVolume() const;
	const T getTotalVolume(const Grid2<bool,T>& isbdry) const;
	const T getVolumeFraction(index_type i, index_type j) const { return MacGridLibrary<T>::heaviside(phi(i,j) * hInv, 0.0f, 1.0f); }

	void getPhi(Grid2<T,T>& _phi, typename MacGrid2D<T>::FIELD_TYPE type) const { getPhi(_phi, phi, type); }
	static void getPhi(Grid2<T,T>&, const Grid2<T,T>&, typename MacGrid2D<T>::FIELD_TYPE type);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawLiquid(const RGBA &color = RGBA::blue()) const;
	void glDrawInterfaceEdges(const RGBA &color = RGBA::blue(), GLfloat width = 3.0f) const;

	void glDrawPhi(const RGBA &colorP = RGBA::red(), const RGBA &colorN = RGBA::white(), T band = 5.0f, GLfloat pointsize = 3.0f) const;	

	static void marchingTriangle(vector<Triangle2D<T>>&, const Grid2<T,T>& phi, const Grid2<bool,T>& isbdry);
	static void marchingInterface(vector<Segment2D<T>>&, const Grid2<T,T>& phi, const Grid2<bool,T>& isbdry);

public : 
	//////////////////////////////////////////////////////////////////////////
	// Variable for distance
	Grid2<T,T>	phi;

	//////////////////////////////////////////////////////////////////////////
	// Variable for Visualization	
	mutable vector<Triangle2D<T>> renderingTriangles;
	mutable vector<Segment2D<T>> renderingSegments;
};

template<typename T>
Levelset2D<T>::Levelset2D() {}

template<typename T>
void Levelset2D<T>::set(index_type nX, index_type nY, T h) {
	MacGrid2D::set(nX, nY, h);

	const Vector2<T> cOffset(T(-0.5), T(-0.5));
	phi.init(nX+2, nY+2, h, h, cOffset); phi.zero();
}

template<typename T>
void Levelset2D<T>::setPhi(const Grid2<T,T>& _phi) {
	phi = _phi; 	
}

template<typename T>
void Levelset2D<T>::setPhi(const DistanceField2D<T>& field) {
	for(index_type i=0;i<phi.nX;++i)	{
		for(index_type j=0;j<phi.nY;++j)	{
			Vector2<T> pos = phi.pos(i,j);
			phi(i,j) = field.getPhi(pos);
		}
	}
	phi.setNeumannBoundary();
}

};		// namespace geometryspace

#include "Levelset2D_Utility.h"
#include "Levelset2D_glDraw.h"

#endif
