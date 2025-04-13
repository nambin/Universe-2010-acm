#pragma once

#include "../MacGrid/MacGrid3D.h"

namespace geometryspace {

template<typename T>
void Levelset3D<T>::createRenderingPrimitive() {
	createRenderingPrimitive(arithmetic::Grid3<bool,T>());
}

template<typename T>
void Levelset3D<T>::createRenderingPrimitive(const arithmetic::Grid3<bool,T>& isbdry) {
	phi.setNeumannBoundary();

	Grid3<T,T> phiNode(nX+1, nY+1, nZ+1, h, h, h, Vector3<T>(0));
	for(index_type i=0;i<phiNode.nX;++i)	
		for(index_type j=0;j<phiNode.nY;++j)	
			for(index_type k=0;k<phiNode.nZ;++k)	
				phiNode(i,j,k) = 0.125f * (phi(i,j,k) + phi(i,j+1,k) + phi(i,j,k+1) + phi(i,j+1,k+1) 
										 + phi(i+1,j,k) + phi(i+1,j+1,k) + phi(i+1,j,k+1) + phi(i+1,j+1,k+1));

	marchingCube(mesh, phiNode, isbdry);
}
//
//template<typename T>
//void Levelset3D<T>::marchingCube(TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>& mesh, const Grid3<T,T>& phi, const Grid3<bool,T>&) {
//
//	const index_type nX = phi.nX-1, nY = phi.nY-1, nZ = phi.nZ-1;
//	const T h = phi.dx, dx = phi.dx, dy = phi.dy, dz = phi.dz;
//
//	const MacGrid3D<T> macgrid(nX,nY,nZ,h);
//
//	Grid3<int,T> newPtOnXEdge, newPtOnYEdge, newPtOnZEdge;
//
//	newPtOnXEdge.init(nX, nY+1, nZ+1, h, h, h);	newPtOnXEdge.set(-1);
//	newPtOnYEdge.init(nX+1, nY, nZ+1, h, h, h);	newPtOnYEdge.set(-1);
//	newPtOnZEdge.init(nX+1, nY+1, nZ, h, h, h);	newPtOnZEdge.set(-1);
//
//	mesh.release();
//
//	for(index_type i=0;i<phi.nX-1;++i)	{
//		for(index_type j=0;j<phi.nY-1;++j)	{
//			for(index_type k=0;k<phi.nZ-1;++k)	{
//
//				int iFlagIndex = 0, iEdgeFlags = 0;
//				index_type nodes[8];
//
//				macgrid.getEightNodesOfCell(i,j,k,nodes);
//
//				for(int iVertex=0;iVertex<8;iVertex++) { if(phi[nodes[iVertex]]<=0) { iFlagIndex |= 1<<iVertex; } } // Which vertices are inside? If i-th vertex is inside, mark '1' at i-th bit of 'iFlagIndex'.
//
//				if(iFlagIndex==0 || iFlagIndex==255) 
//					continue; 
//
//				int *newPtOnEdge[12];
//
//				macgrid.getEdgesInfoOfCell(i,j,k, newPtOnEdge, newPtOnXEdge, newPtOnYEdge, newPtOnZEdge);
//				iEdgeFlags = CubeEdgeFlags3D[iFlagIndex]; // Which edges intersect the surface? If i-th edge intersects the surface, mark '1' at i-th bit of 'iEdgeFlags' 
//
//				Vector3<T> pt[12];
//
//				const T EPS = 1.e-3;
//				for( int iEdge=0; iEdge<12; iEdge++ ) { // Find the point of intersection of the surface per each edge
//
//					if( iEdgeFlags & (1<<iEdge) ) { // If there is an intersection on this edge
//
//						int e0 = EdgeConnection3D[iEdge][0];
//						int e1 = EdgeConnection3D[iEdge][1];
//
//						T alpha = abs( phi[nodes[e0]] ) / abs( phi[nodes[e1]] - phi[nodes[e0]] );
//
//						if(alpha < EPS)				alpha = EPS;
//						else if(alpha > 1.0f - EPS)	alpha = 1.0f - EPS;
//
//						pt[iEdge][0] = dx * ( i + VertexOffset3D[e0][0] + alpha * EdgeDirection3D[iEdge][0] );
//						pt[iEdge][1] = dy * ( j + VertexOffset3D[e0][1] + alpha * EdgeDirection3D[iEdge][1] );
//						pt[iEdge][2] = dz * ( k + VertexOffset3D[e0][2] + alpha * EdgeDirection3D[iEdge][2] );
//
//						if( *newPtOnEdge[iEdge] == -1 ) {
//							mesh.posVerts.push_back( pt[iEdge] + phi.origin );
//							*newPtOnEdge[iEdge] = mesh.numVertices++;
//						}
//					}
//				}
//
//				int threeVerts[3];
//
//				for( int iTriangle=0; iTriangle<5; ++iTriangle ) {
//
//					// If there isn't any triangle to be drawn, escape this loop.
//					if( TriangleConnectionTableMC3D[iFlagIndex][3*iTriangle] < 0 ) break;
//
//					for( int iCorner=0; iCorner<3; iCorner++ ) {
//
//						int iVertex = TriangleConnectionTableMC3D[iFlagIndex][3*iTriangle+iCorner];						
//						threeVerts[iCorner] = *newPtOnEdge[iVertex];
//
//					}
//
//					if(threeVerts[0] == -1 || threeVerts[1] == -1 || threeVerts[2] == -1)	{
//						std::cout << "Error in MC reconstruction !" << std::endl;
//						std::cout << i << " " << j << " " << k << " " << iFlagIndex << std::endl;
//					}
//
//					mesh.v2f.push_back( threeVerts[0] );
//					mesh.v2f.push_back( threeVerts[1] );
//					mesh.v2f.push_back( threeVerts[2] );
//					mesh.numFaces++;
//				}
//			}
//		}
//	}
//
//	mesh.makeFaceNormals();
//	mesh.makeVertexNormals();
//}

template<typename T>
void Levelset3D<T>::marchingCube(TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>& mesh, const Grid3<T,T>& phi, const Grid3<bool,T>& isbdry) {

	vector<Triangle3D<T>> tmp;
	Vector3<T> org, h(phi.dx, phi.dy, phi.dz);

	mesh.release();

	for(index_type i=1;i<phi.nX;++i)	{
		for(index_type j=1;j<phi.nY;++j)	{
			for(index_type k=1;k<phi.nZ;++k)	{

				if(isbdry.is_empty() == false && isbdry(i,j,k))
					continue;

				T _phi[8];

				_phi[0] = phi(i-1,j-1,k-1);	_phi[1] = phi(i  ,j-1,k-1);	_phi[2] = phi(i  ,j-1,k  );	_phi[3] = phi(i-1,j-1,k  );
				_phi[4] = phi(i-1,j  ,k-1);	_phi[5] = phi(i  ,j  ,k-1);	_phi[6] = phi(i  ,j  ,k  );	_phi[7] = phi(i-1,j  ,k  );				

				org = phi.pos(i-1,j-1,k-1);			

				if(isbdry.is_empty() == false)
					contouring(tmp, _phi, org, h, isbdry(i-1,j,k), isbdry(i+1,j,k)
											    , isbdry(i,j-1,k), isbdry(i,j+1,k)
											    , isbdry(i,j,k-1), isbdry(i,j,k+1));
				else 
					contouring(tmp, _phi, org, h, false, false, false, false, false, false);

				for(size_t m=0;m<tmp.size();++m) {

					mesh.posVerts.push_back(tmp[m].v[0]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1));
					mesh.posVerts.push_back(tmp[m].v[1]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1));
					mesh.posVerts.push_back(tmp[m].v[2]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1));

				}				
			}
		}
	}

	mesh.numVertices = unsigned int(mesh.posVerts.size());
	mesh.numFaces = unsigned int(mesh.v2f.size() / 3);

	Grid3<Vector3<T>,T> normals; normals.resize(phi);
	MacGridLibrary<T>::getNormals(normals, phi);

	mesh.normals.resize(mesh.posVerts.size());

#ifdef USE_OPENMP	
	int my_end = int(mesh.posVerts.size());
	int i;

#pragma omp parallel for
#else	
	size_t my_end = mesh.posVerts.size();
	size_t i;
#endif
	for(i=0;i<my_end;++i) 
		mesh.normals[i] = unitVector(normals(mesh.posVerts[i]));
}

template<typename T>
const T Levelset3D<T>::getTotalVolume() const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			for(index_type k=1;k<phi.nZ-1;++k)
				vol += getVolumeFraction(i,j,k);
	return vol * h * h * h;
}

template<typename T>
const T Levelset3D<T>::getTotalVolume(const Grid3<bool,T>& isbdry) const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			for(index_type k=1;k<phi.nZ-1;++k)
				if(isbdry(i,j,k) == false)
					vol += getVolumeFraction(i,j,k);
	return vol * h * h * h;
}

template<typename T>
void Levelset3D<T>::getPhi(Grid3<T,T>& _phi, const Grid3<T,T>& phi, typename MacGrid3D<T>::FIELD_TYPE _type) {
	if(_type == MacGrid3D<T>::TYPE_U) {
		assert(_phi.nX == phi.nX-1 && _phi.nY == phi.nY && _phi.nZ == phi.nZ);

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
				for(index_type k=0;k<_phi.nZ;++k) 
					_phi(i,j,k) = T(0.5) * (phi(i,j,k) + phi(i+1,j,k));
	}
	else if(_type == MacGrid3D<T>::TYPE_V) {
		assert(_phi.nX == phi.nX && _phi.nY == phi.nY-1 && _phi.nZ == phi.nZ);

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
				for(index_type k=0;k<_phi.nZ;++k) 
					_phi(i,j,k) = T(0.5) * (phi(i,j,k) + phi(i,j+1,k));
	}
	else if(_type == MacGrid3D<T>::TYPE_W) {
		assert(_phi.nX == phi.nX && _phi.nY == phi.nY && _phi.nZ == phi.nZ-1);

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
				for(index_type k=0;k<_phi.nZ;++k) 
					_phi(i,j,k) = T(0.5) * (phi(i,j,k) + phi(i,j,k+1));
	}
	else if(_type == MacGrid3D<T>::TYPE_S) {
		assert(_phi.nX == phi.nX && _phi.nY == phi.nY && _phi.nZ == phi.nZ);

		_phi = phi;
	}
}

};			// geometryspace