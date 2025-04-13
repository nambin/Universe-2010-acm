#ifndef _GEOMETRYSPACE_MACGRID_3D_H_
#define _GEOMETRYSPACE_MACGRID_3D_H_

#include "../../Arithmetic/MathLibrary/Interpolation.h"
#include "../../VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h"
#include "../../VirtualSpace/DrawLibrary/glDrawLibrary.h"
#include "../../Arithmetic/GridLibrary/Grid3.h"

namespace geometryspace {

using basis::RGBA;
using arithmetic::Vector3;

template<typename T>
class MacGrid3D	{
public : 
	typedef typename arithmetic::Grid3<T,T>::index_type index_type;

	enum FIELD_TYPE { TYPE_U = 0, TYPE_V, TYPE_W, TYPE_S, TYPE_UNDEFINED };

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	MacGrid3D()	{}
	MacGrid3D(index_type nX, index_type nY, index_type nZ, T h) { set(nX, nY, nZ, h); }
	virtual ~MacGrid3D() {}

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	virtual void set(index_type nX, index_type nY, index_type nZ, T h);

	const bool importBinary(std::ifstream&);
	const bool exportBinary(std::ofstream&) const;	

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawBoundingBox( const RGBA &color = RGBA::white(), GLfloat width = 1.5f) const;

public : 
	T				Lx, Ly, Lz;
	index_type		nX, nY, nZ;
	index_type		numU, numV, numW, numS;
	index_type		numNodes, numCells;
	index_type		numXEdges, numYEdges, numZEdges;

	T				dx, dy, dz;
	T				dxInv, dyInv, dzInv;
	T				dx2, dy2, dz2;
	T				h, hInv;
	T				h2;

public :
	inline index_type iU(index_type i, index_type j, index_type k) const { return k + (nZ+2) * (j + (nY+2) * i); }
	inline index_type iV(index_type i, index_type j, index_type k) const { return k + (nZ+2) * (j + (nY+1) * i); }
	inline index_type iW(index_type i, index_type j, index_type k) const { return k + (nZ+1) * (j + (nY+2) * i); }
	inline index_type iS(index_type i, index_type j, index_type k) const { return k + (nZ+2) * (j + (nY+2) * i); }
	inline index_type iNode(index_type i, index_type j, index_type k) const { return k + (nZ+1) * (j + (nY+1) * i);	}
	inline index_type iXEdge(index_type i, index_type j, index_type k) const { return k + (nZ+1) * (j + (nY+1) * i);}
	inline index_type iYEdge(index_type i, index_type j, index_type k) const { return k + (nZ+1) * (j + (nY) * i);	}
	inline index_type iZEdge(index_type i, index_type j, index_type k) const { return k + (nZ) * (j + (nY+1) * i);	}

	inline index_type iS2U(index_type i, index_type j, index_type k, index_type dir) const { return iU(i+dir-1,j,k);	}
	inline index_type iS2V(index_type i, index_type j, index_type k, index_type dir) const { return iV(i,j+dir-1,k);	}
	inline index_type iS2W(index_type i, index_type j, index_type k, index_type dir) const { return iW(i,j,k+dir-1);	}
	inline index_type iU2S(index_type i, index_type j, index_type k, index_type dir) const { return iS(i+dir,j,k);		}
	inline index_type iV2S(index_type i, index_type j, index_type k, index_type dir) const { return iS(i,j+dir,k);		}
	inline index_type iW2S(index_type i, index_type j, index_type k, index_type dir) const { return iS(i,j,k+dir);		}

	Vector3<T> posU(index_type i, index_type j, index_type k) const { return Vector3<T>(i*dx, (j-0.5f)*dy, (k-0.5f)*dz);		}
	Vector3<T> posV(index_type i, index_type j, index_type k) const { return Vector3<T>((i-0.5f)*dx, j*dy, (k-0.5f)*dz);		}
	Vector3<T> posW(index_type i, index_type j, index_type k) const { return Vector3<T>((i-0.5f)*dx, (j-0.5f)*dy, k*dz);		}
	Vector3<T> posS(index_type i, index_type j, index_type k) const { return Vector3<T>((i-0.5f)*dx, (j-0.5f)*dy, (k-0.5f)*dz);	}
	Vector3<T> posNode(index_type i, index_type j, index_type k) const { return Vector3<T>(i*dx, j*dy, k*dz);		}
	Vector3<T> posNode(index_type iNode) const { 
		index_type i = iNode / ((nY+1)*(nZ+1));
		index_type j = (iNode - ((nY+1)*(nZ+1)) * i) / (nZ+1);
		index_type k = iNode % (nZ+1);		

		return Vector3<T>(i*dx, j*dy, k*dz);		
	}

	const bool isInsideDomain(const Vector3<T> &pos) const {
		if(pos[0] < 0 || pos[0] > Lx 
		|| pos[1] < 0 || pos[1] > Ly
		|| pos[2] < 0 || pos[2] > Lz) 	
			return false;
		else	
			return true;
	}

	const geometryspace::AABB<T> getAABB() const {
		return geometryspace::AABB<T>(Vector3<T>(0), Vector3<T>(Lx, Ly, Lz));
	}

	void constrainedBdry( Vector3<T> &pos ) const {
		static const T eps = 1.0e-8;
		if ( pos[0] < eps ) pos[0] = eps; else if ( pos[0] > Lx-eps ) pos[0] = Lx - eps;
		if ( pos[1] < eps ) pos[1] = eps; else if ( pos[1] > Ly-eps ) pos[1] = Ly - eps;
		if ( pos[2] < eps ) pos[2] = eps; else if ( pos[2] > Lz-eps ) pos[2] = Lz - eps;
	}

	void getIdxU( const Vector3<T>& pos, index_type idxU[3]) const {
		idxU[0] = (index_type)(pos[0]*dxInv); 
		idxU[1] = (index_type)(pos[1]*dyInv+0.5f);
		idxU[2] = (index_type)(pos[2]*dzInv+0.5f);
		if ( idxU[0] == nX   ) idxU[0]--;	
		if ( idxU[1] == nY+1 ) idxU[1]--;	
		if ( idxU[2] == nZ+1 ) idxU[2]--;
	}
	void getIdxV( const Vector3<T>& pos, index_type idxV[3]) const {
		idxV[0] = (index_type)(pos[0]*dxInv+0.5f); 
		idxV[1] = (index_type)(pos[1]*dyInv);
		idxV[2] = (index_type)(pos[2]*dzInv+0.5f);
		if ( idxV[0] == nX+1 ) idxV[0]--;	
		if ( idxV[1] == nY   ) idxV[1]--;	
		if ( idxV[2] == nZ+1 ) idxV[2]--;
	}
	void getIdxW( const Vector3<T> &pos, index_type idxW[3] ) const {
		idxW[0] = (index_type)(pos[0]*dxInv+0.5f); 
		idxW[1] = (index_type)(pos[1]*dyInv+0.5f);
		idxW[2] = (index_type)(pos[2]*dzInv);
		if ( idxW[0] == nX+1 ) idxW[0]--; 
		if ( idxW[1] == nY+1 ) idxW[1]--;	
		if ( idxW[2] == nZ   ) idxW[2]--;
	}
	void getIdxS( const Vector3<T> &pos, index_type idxS[3] ) const {
		idxS[0] = (index_type)(pos[0]*dxInv+0.5f); 
		idxS[1] = (index_type)(pos[1]*dyInv+0.5f);
		idxS[2] = (index_type)(pos[2]*dzInv+0.5f);
		if ( idxS[0] == nX+1 ) idxS[0]--; 
		if ( idxS[1] == nY+1 ) idxS[1]--;	
		if ( idxS[2] == nZ+1 ) idxS[2]--;
	}
	void getIdxNode( const Vector3<T> &pos, index_type idxNode[3] ) const {
		idxNode[0] = (index_type)(pos[0]*dxInv);
		idxNode[1] = (index_type)(pos[1]*dyInv);
		idxNode[2] = (index_type)(pos[2]*dzInv);
		if ( idxNode[0] == nX ) idxNode[0]--; 
		if ( idxNode[1] == nY ) idxNode[1]--;	
		if ( idxNode[2] == nZ ) idxNode[2]--;
	}
	index_type idxU( const Vector3<T> &pos ) const {
		index_type idx[3];
		getIdxU(pos, idx);
		return iU(idx[0],idx[1],idx[2]);
	}
	index_type idxV( const Vector3<T> &pos ) const {
		index_type idx[3];
		getIdxV(pos, idx);
		return iV(idx[0],idx[1],idx[2]);
	}
	index_type idxW( const Vector3<T> &pos ) const {
		index_type idx[3];
		getIdxW(pos, idx);
		return iW(idx[0],idx[1],idx[2]);
	}
	index_type idxS( const Vector3<T> &pos ) const {
		index_type idx[3];
		getIdxS(pos, idx);
		return iS(idx[0],idx[1],idx[2]);
	}
	index_type idxNode( const Vector3<T> &pos ) const {
		index_type idx[3];
		getIdxNode(pos, idx);
		return iNode(idx[0],idx[1],idx[2]);
	}

	void getTwelveEdgesOfCell( index_type i, index_type j, index_type k, index_type xEdge[4], index_type yEdge[4], index_type zEdge[4] ) const {
		xEdge[0] = iXEdge(i  ,j  ,k  );	// edge 0(0-1)
		zEdge[0] = iZEdge(i+1,j  ,k  );	// edge 1(1-2)
		xEdge[1] = iXEdge(i  ,j  ,k+1);	// edge 2(2-3)
		zEdge[1] = iZEdge(i  ,j  ,k  );	// edge 3(3-0)

		xEdge[2] = iXEdge(i  ,j+1,k  );	// edge 4(4-5)
		zEdge[2] = iZEdge(i+1,j+1,k  );	// edge 5(5-6)
		xEdge[3] = iXEdge(i  ,j+1,k+1);	// edge 6(6-7)
		zEdge[3] = iZEdge(i  ,j+1,k  );	// edge 7(7-4)

		yEdge[0] = iYEdge(i  ,j  ,k  );	// edge 8(0-4)
		yEdge[1] = iYEdge(i+1,j  ,k  );	// edge 9(1-5)
		yEdge[2] = iYEdge(i+1,j  ,k+1);	// edge10(2-6)
		yEdge[3] = iYEdge(i  ,j  ,k+1);	// edge11(3-7)
	}

	void getEightNodesOfCell( index_type i, index_type j, index_type k, index_type nodes[8] ) const {
		nodes[0] = iNode(i  ,j  ,k  );	// node0
		nodes[1] = iNode(i+1,j  ,k  );	// node1
		nodes[2] = iNode(i+1,j  ,k+1);	// node2
		nodes[3] = iNode(i  ,j  ,k+1);	// node3
		nodes[4] = iNode(i  ,j+1,k  );	// node4
		nodes[5] = iNode(i+1,j+1,k  );	// node5
		nodes[6] = iNode(i+1,j+1,k+1);	// node6
		nodes[7] = iNode(i  ,j+1,k+1);	// node7
	}

	void getEdgeInfo(index_type idx, int *&newPtOnEdge, Grid3<int,T> &newPtOnWEdge)	const {
		newPtOnEdge = &(newPtOnWEdge.data[idx]);
	}
	void getEdgesInfoOfCell( index_type i, index_type j, index_type k, int *newPtOnEdge[12]
						   , Grid3<int,T>& newPtOnXEdge, Grid3<int,T>& newPtOnYEdge, Grid3<int,T>& newPtOnZEdge ) const {

		index_type xEdge[4], yEdge[4], zEdge[4];
		getTwelveEdgesOfCell( i, j, k, xEdge, yEdge, zEdge );

		for(std::size_t i=0;i<12;i++)	{
			if( i > 7 )				getEdgeInfo( yEdge[i-8], newPtOnEdge[i], newPtOnYEdge );
			else if ( i % 2 == 0 ) 	getEdgeInfo( xEdge[(size_t)(i/2)], newPtOnEdge[i], newPtOnXEdge );
			else					getEdgeInfo( zEdge[(size_t)(i/2)], newPtOnEdge[i], newPtOnZEdge );
		}
	}

	inline T getS(const Vector3<T> &_pos, const T *_s) const {

		index_type r = (index_type)(_pos[0]*dxInv+0.5f); 
		index_type s = (index_type)(_pos[1]*dyInv+0.5f);
		index_type t = (index_type)(_pos[2]*dzInv+0.5f);

		if ( r == nX+1 ) r--;	
		if ( s == nY+1 ) s--;	
		if ( t == nZ+1 ) t--;

		Vector3<T> corner = posS(r,s,t);
		T alpha = ( _pos[0] - corner[0] ) * dxInv;
		T beta  = ( _pos[1] - corner[1] ) * dyInv;
		T gamma = ( _pos[2] - corner[2] ) * dzInv;

		return (1.f-alpha)*( (1.f-beta)	* ((1.f-gamma)*_s[iS(r,s,t)]		+ gamma*_s[iS(r,s,t+1)])
								+ beta	* ((1.f-gamma)*_s[iS(r,s+1,t)]		+ gamma*_s[iS(r,s+1,t+1)])   )
				   + alpha*( (1.f-beta) * ((1.f-gamma)*_s[iS(r+1,s,t)]		+ gamma*_s[iS(r+1,s,t+1)]) 
								+ beta	* ((1.f-gamma)*_s[iS(r+1,s+1,t)]	+ gamma*_s[iS(r+1,s+1,t+1)]) );
	}
	inline T getU(const Vector3<T> &_pos, const T *_u) const {

		index_type r = (index_type)(_pos[0]*dxInv); 
		index_type s = (index_type)(_pos[1]*dyInv+0.5f);
		index_type t = (index_type)(_pos[2]*dzInv+0.5f);

		if ( r == nX ) r--;	if ( s == nY+1 ) s--;	if ( t == nZ+1 ) t--;

		Vector3<T> corner = posU(r,s,t);
		T alpha = ( _pos[0] - corner[0] ) * dxInv;
		T beta  = ( _pos[1] - corner[1] ) * dyInv;
		T gamma = ( _pos[2] - corner[2] ) * dzInv;

		return (1.f-alpha)*( (1.f-beta)	* ((1.f-gamma)*_u[iU(r,s,t)]		+ gamma*_u[iU(r,s,t+1)])
								+ beta	* ((1.f-gamma)*_u[iU(r,s+1,t)]		+ gamma*_u[iU(r,s+1,t+1)])   )
				   + alpha*( (1.f-beta) * ((1.f-gamma)*_u[iU(r+1,s,t)]		+ gamma*_u[iU(r+1,s,t+1)]) 
								+ beta	* ((1.f-gamma)*_u[iU(r+1,s+1,t)]	+ gamma*_u[iU(r+1,s+1,t+1)]) );
	}
	inline T getV(const Vector3<T> &_pos, const T *_v) const {

		index_type r = (index_type)(_pos[0]*dxInv+0.5f); 
		index_type s = (index_type)(_pos[1]*dyInv);
		index_type t = (index_type)(_pos[2]*dzInv+0.5f);

		if ( r == nX+1 ) r--;	if ( s == nY ) s--;	if ( t == nZ+1 ) t--;

		Vector3<T> corner = posV(r,s,t);
		T alpha = ( _pos[0] - corner[0] ) * dxInv;
		T beta  = ( _pos[1] - corner[1] ) * dyInv;
		T gamma = ( _pos[2] - corner[2] ) * dzInv;

		return (1.f-alpha)*( (1.f-beta)	* ((1.f-gamma)*_v[iV(r,s,t)]		+ gamma*_v[iV(r,s,t+1)])
								+ beta	* ((1.f-gamma)*_v[iV(r,s+1,t)]		+ gamma*_v[iV(r,s+1,t+1)])   )
				   + alpha*( (1.f-beta) * ((1.f-gamma)*_v[iV(r+1,s,t)]		+ gamma*_v[iV(r+1,s,t+1)]) 
								+ beta	* ((1.f-gamma)*_v[iV(r+1,s+1,t)]	+ gamma*_v[iV(r+1,s+1,t+1)]) );

	}
	inline T getW(const Vector3<T> &_pos, const T *_w) const {

		index_type r = (index_type)(_pos[0]*dxInv+0.5f); 
		index_type s = (index_type)(_pos[1]*dyInv+0.5f);
		index_type t = (index_type)(_pos[2]*dzInv);

		if ( r == nX+1 ) r--;	if ( s == nY+1 ) s--;	if ( t == nZ ) t--;

		Vector3<T> corner = posW(r,s,t);
		T alpha = ( _pos[0] - corner[0] ) * dxInv;
		T beta  = ( _pos[1] - corner[1] ) * dyInv;
		T gamma = ( _pos[2] - corner[2] ) * dzInv;

		return (1.f-alpha)*( (1.f-beta)	* ((1.f-gamma)*_w[iW(r,s,t)]		+ gamma*_w[iW(r,s,t+1)])
								+ beta	* ((1.f-gamma)*_w[iW(r,s+1,t)]		+ gamma*_w[iW(r,s+1,t+1)])   )
				   + alpha*( (1.f-beta) * ((1.f-gamma)*_w[iW(r+1,s,t)]		+ gamma*_w[iW(r+1,s,t+1)]) 
								+ beta	* ((1.f-gamma)*_w[iW(r+1,s+1,t)]	+ gamma*_w[iW(r+1,s+1,t+1)]) );
	}
};

};		// geometryspace

#include "MacGrid3D_Functions.h"

#endif