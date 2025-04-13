#pragma once

#include "../../../Basis/vectorLibrary/vectorLibrary.h"

namespace geometryspace {

template<typename T>
void Levelset2D<T>::createRenderingPrimitive() {
	createRenderingPrimitive(arithmetic::Grid2<bool,T>());
}

template<typename T>
void Levelset2D<T>::createRenderingPrimitive(const arithmetic::Grid2<bool,T>& isbdry) {
	Grid2<T,T> phiNode(nX+1, nY+1, h, h, Vector2<T>(0));

	phi.setNeumannBoundary();
	for(index_type i=0;i<phiNode.nX;++i) 
		for(index_type j=0;j<phiNode.nY;++j) 
			phiNode(i,j) = 0.25f * (phi(i,j) + phi(i+1,j) + phi(i,j+1) + phi(i+1,j+1));

	marchingTriangle(renderingTriangles, phiNode, isbdry);
	marchingInterface(renderingSegments, phiNode, isbdry);
}

template<typename T>
void Levelset2D<T>::marchingInterface(vector<Segment2D<T>>& segments, const Grid2<T,T>& phiNode, const Grid2<bool,T>& isbdry) {
	segments.clear();

	vector<Segment2D<T>> tmp;
	Vector2<T> org, h(phiNode.dx, phiNode.dy);

	for(index_type i=1;i<phiNode.nX;++i)	{
		for(index_type j=1;j<phiNode.nY;++j)	{

			if(isbdry.is_empty() == false && isbdry(i,j))
				continue;

			T phi[4];

			phi[0] = phiNode(i-1, j-1);	phi[1] = phiNode(i  , j-1);
			phi[2] = phiNode(i  , j  );	phi[3] = phiNode(i-1, j  );

			org = phiNode.pos(i-1, j-1);	

			if(isbdry.is_empty() == false)
				contouring(tmp, phi, org, h, isbdry(i-1,j), isbdry(i+1,j)
										   , isbdry(i,j-1), isbdry(i,j+1));
			else 
				contouring(tmp, phi, org, h, false, false, false, false);

			basis::vectorLibrary::append(segments, tmp);
		}
	}	
}

template<typename T>
void Levelset2D<T>::marchingTriangle(vector<Triangle2D<T>>& renderingTriangles, const Grid2<T,T> &phiNode, const Grid2<bool,T>& isbdry) {
	renderingTriangles.clear();

	vector<Triangle2D<T>> tmp;
	Vector2<T> org, h(phiNode.dx, phiNode.dy);

	for(index_type i=1;i<phiNode.nX;++i)	{
		for(index_type j=1;j<phiNode.nY;++j)	{

			if(isbdry.is_empty() == false && isbdry(i,j))
				continue;

			T phi[4];

			phi[0] = phiNode(i-1, j-1);	phi[1] = phiNode(i  , j-1);
			phi[2] = phiNode(i  , j  );	phi[3] = phiNode(i-1, j  );

			org = phiNode.pos(i-1, j-1);			
			
			contouring(tmp, phi, org, h);

			basis::vectorLibrary::append(renderingTriangles, tmp);
		}
	}	
}

template<typename T>
const T Levelset2D<T>::getTotalVolume() const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			vol += getVolumeFraction(i,j);
	return vol * h * h;
}

template<typename T>
const T Levelset2D<T>::getTotalVolume(const Grid2<bool,T>& isbdry) const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			if(isbdry(i,j) == false)
				vol += getVolumeFraction(i,j);
	return vol * h * h;
}

template<typename T>
void Levelset2D<T>::getPhi(Grid2<T,T>& _phi, const Grid2<T,T>& phi, typename MacGrid2D<T>::FIELD_TYPE _type) {
	if(_type == MacGrid2D<T>::TYPE_U) {
		assert(_phi.nX == phi.nX-1 && _phi.nY == phi.nY);

#ifdef USE_OPENMP
		int my_begin = int(0);
		int my_end = int(_phi.nX);
		int i;

#pragma omp parallel for
#else
		index_type my_begin = 0;
		index_type my_end = _phi.nX;
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 		
			for(index_type j=0;j<_phi.nY;++j) 
				_phi(i,j) = T(0.5) * (phi(i,j) + phi(i+1,j));
	}
	else if(_type == MacGrid2D<T>::TYPE_V) {
		assert(_phi.nX == phi.nX && _phi.nY == phi.nY-1);

#ifdef USE_OPENMP
		int my_begin = int(0);
		int my_end = int(_phi.nX);
		int i;

#pragma omp parallel for
#else
		index_type my_begin = 0;
		index_type my_end = _phi.nX;
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 	
			for(index_type j=0;j<_phi.nY;++j) 
				_phi(i,j) = T(0.5) * (phi(i,j) + phi(i,j+1));
	}
	else if(_type == MacGrid2D<T>::TYPE_S) {
		assert(_phi.nX == phi.nX && _phi.nY == phi.nY);

		_phi = phi;
	}
}

};		// geometryspace