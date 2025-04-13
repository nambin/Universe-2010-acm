#ifndef _GEOMETRYSPACE_MACGRID_2D_H_
#define _GEOMETRYSPACE_MACGRID_2D_H_

#include "../../../Basis/FileIO/FileIO.h"
#include "../../VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h"
#include "../../Arithmetic/MathLibrary/Interpolation.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/GridLibrary/Grid2.h"

namespace geometryspace {

using basis::RGBA;
using arithmetic::Vector2;

template<typename T>
class MacGrid2D	{
public : 
	typedef typename arithmetic::Grid2<T,T>::index_type index_type;

	enum FIELD_TYPE { TYPE_U = 0, TYPE_V, TYPE_S, TYPE_UNDEFINED };

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	MacGrid2D()	{}
	MacGrid2D(index_type nX, index_type nY, T h) { set(nX, nY, h); }
	virtual ~MacGrid2D() {}

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	virtual void set(index_type nX, index_type nY, T h);	

	const bool importBinary(std::ifstream&);
	const bool exportBinary(std::ofstream&) const;	

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawBoundingBox( const RGBA &color, GLfloat width ) const;
	void glDrawUniformGrid( const RGBA &color, GLfloat width ) const;	

	void glDrawNode(index_type i, index_type j, const RGBA &color, GLfloat size) const;
	void glDrawNode(index_type i, const RGBA &color, GLfloat size) const;	
	void glDrawEdge(index_type i, const RGBA &color, GLfloat width) const;
	void glDrawCell(index_type i, index_type j, const RGBA &color) const;
	void glDrawCell(index_type i, const RGBA &color) const;	

public : 
	T				Lx, Ly;
	index_type		nX, nY;
	index_type		numU, numV, numS;
	index_type		numNodes, numEdges;

	T				h, hInv, h2;

public :	
	inline const index_type iU(index_type i, index_type j) const		{ return i*(nY+2) + j;			}
	inline const index_type iV(index_type i, index_type j) const		{ return i*(nY+1) + j;			}
	inline const index_type iS(index_type i, index_type j) const		{ return i*(nY+2) + j;			}
	inline const index_type iNode(index_type i, index_type j) const		{ return i*(nY+1) + j;			}
	inline const index_type iXEdge(index_type i, index_type j) const	{ return j*nX + i;				}
	inline const index_type iYEdge(index_type i, index_type j) const	{ return i*nY + j + nX * (nY+1);}

	inline const index_type iS2U( index_type i, index_type j, index_type dir ) const { return iU(i+dir-1,j);}
	inline const index_type iS2V( index_type i, index_type j, index_type dir ) const { return iV(i,j+dir-1);}
	inline const index_type iU2S( index_type i, index_type j, index_type dir ) const { return iS(i+dir,j);	}
	inline const index_type iV2S( index_type i, index_type j, index_type dir ) const { return iS(i,j+dir);	}

	const bool isOuterCell( index_type cell ) const {
		index_type i = cell / (nY+2);
		index_type j = cell % (nY+2);

		if(i==0 || i==nX+1 || j==0 || j==nY+1)	return true;
		else									return false;
	}

	inline const Vector2<T> posU(index_type i, index_type j) const { return Vector2<T>(i*h      , (j-0.5f)*h);	}
	inline const Vector2<T> posV(index_type i, index_type j) const { return Vector2<T>((i-0.5f)*h, j*h	   );	}
	inline const Vector2<T> posS(index_type i, index_type j) const { return Vector2<T>((i-0.5f)*h, (j-0.5f)*h);	}
	inline const Vector2<T> posNode(index_type i, index_type j) const { return Vector2<T>(i*h, j*h);			}
	inline const Vector2<T> posNode(index_type iNode) const { 
		index_type i = iNode / (nY+1);
		index_type j = iNode % (nY+1);
		return Vector2<T>(i*h, j*h);		
	}

	const bool isInsideDomain(const Vector2<T> &pos) const {
		if(pos[0] < 0 || pos[0] > Lx 
		|| pos[1] < 0 || pos[1] > Ly) 	
			return false;
		else
			return true;
	}

	const geometryspace::AABB2D<T> getAABB() const {
		return geometryspace::AABB2D<T>(Vector2<T>(0), Vector2<T>(Lx, Ly));
	}

	void constrainedBdry( Vector2<T> &pos ) const {
		static const T eps = 1.0e-8;
		if ( pos[0] < eps ) pos[0] = eps; else if ( pos[0] > Lx-eps ) pos[0] = Lx - eps;
		if ( pos[1] < eps ) pos[1] = eps; else if ( pos[1] > Ly-eps ) pos[1] = Ly - eps;
	}

	inline index_type idxU( const Vector2<T> &pos ) const {
		index_type r = (index_type)(pos[0]*hInv); 
		index_type s = (index_type)(pos[1]*hInv+0.5f);
		if ( r == nX ) r--;		if ( s == nY+1 ) s--;
		return iU(r,s);
	}
	inline index_type idxV( const Vector2<T> &pos ) const {
		index_type r = (index_type)(pos[0]*hInv+0.5f); 
		index_type s = (index_type)(pos[1]*hInv);
		if ( r == nX+1 ) r--;	if ( s == nY ) s--;
		return iV(r,s);
	}
	inline index_type idxS( const Vector2<T> &pos ) const {
		index_type r = (index_type)(pos[0]*hInv+0.5f); 
		index_type s = (index_type)(pos[1]*hInv+0.5f);
		if ( r == nX+1 ) r--;	if ( s == nY+1 ) s--;
		return iS(r,s);
	}
	inline index_type idxNode( const Vector2<T> &pos ) const {
		index_type r = (index_type)(pos[0]*hInv);
		index_type s = (index_type)(pos[1]*hInv);
		if ( r == nX ) r--;		if ( s == nY ) s--;
		return iNode(r,s);
	}

	inline index_type idxN2N( index_type i, index_type j, index_type dir ) const {
		switch(dir)	{
		case 0 : if(i==nX)	return -1;
				 else		return iNode(i+1,j);
		case 1 : if(j==nY)	return -1;
				 else		return iNode(i,j+1);
		case 2 : if(i==0)	return -1;
				 else		return iNode(i-1,j);
		case 3 : if(j==0)	return -1;
				 else		return iNode(i,j-1);
		default : return -1;
		}		
	}
	inline index_type idxN2N( index_type iNode, index_type dir ) const {
		index_type i = iNode / (nY+1);
		index_type j = iNode % (nY+1);
		return idxN2N(i,j,dir);
	}
	inline index_type idxE2N( index_type i, index_type j, index_type dir ) const {	
		switch(dir)	{
		case 0 : if(i==nX)	return -1;
				 else		return nX * j + i;
		case 1 : if(j==nY)	return -1;
				 else		return i * nY + j + nX * (nY+1);
		case 2 : if(i==0)	return -1;
				 else		return nX * j + i-1;
		case 3 : if(j==0)	return -1;
				 else		return i * nY + j-1 + nX * (nY+1);
		default : return -1;
		}		
	}
	inline index_type idxE2N( index_type iNode, index_type dir ) const {
		index_type i = iNode / (nY+1);
		index_type j = iNode % (nY+1);
		return idxE2N(i,j,dir);
	}
	inline index_type idxNodeEdge( index_type iEdge, index_type dir ) const {
		index_type i,j;
		if(iEdge < nX * (nY+1))	{
			i = iEdge % nX;
			j = iEdge / nX;
			if(dir == 0)	return iNode(i, j);
			if(dir == 1)	return iNode(i+1, j);
		}
		else	{
			i = (iEdge - nX * (nY+1)) / nY;
			j = (iEdge - nX * (nY+1)) % nY;
			if(dir == 0)	return iNode(i, j);
			if(dir == 1)	return iNode(i, j+1);
		}		
		return -1;
	}
	inline index_type idxS2E( index_type iEdge, index_type dir ) const {
		index_type i,j;
		if(iEdge < nX * (nY+1))	{
			i = iEdge % nX;
			j = iEdge / nX;
			if(dir == 0)	return iS(i+1, j);
			if(dir == 1)	return iS(i+1, j+1);
			return -1;
		}
		else	{
			i = (iEdge - nX * (nY+1)) / nY;
			j = (iEdge - nX * (nY+1)) % nY;
			if(dir == 0)	return iS(i, j+1);
			if(dir == 1)	return iS(i+1, j+1);
			return -1;
		}
	}
	inline index_type idxS2N( index_type i, index_type j, index_type dir) const {
		switch(dir)	{
		case 0 : return iS(i  ,j  );
		case 1 : return iS(i+1,j  );
		case 2 : return iS(i+1,j+1);
		case 3 : return iS(i  ,j+1);
		default : return -1;
		}		
	}
	inline index_type idxS2N( index_type iNode, index_type dir ) const {
		index_type i = iNode / (nY+1);
		index_type j = iNode % (nY+1);
		return idxS2N(i,j,dir);
	}
	inline index_type idxEdgeCell( index_type i, index_type j, index_type dir ) const {
		if((i == 0 && dir != 1) || (i == nX+1 && dir != 3))	return -1;
		if((j == 0 && dir != 2) || (j == nY+1 && dir != 0))	return -1;

		switch(dir) {
		case 0 : return ((j-1)*nX + i-1);
		case 2 : return (j*nX + i-1);
		case 3 : return ((i-1)*nY + j-1 + nX * (nY+1));
		case 1 : return (i*nY + j-1 + nX * (nY+1));
		default : return -1;
		}		
	}
	inline index_type idxEdgeCell( index_type iSell, index_type dir) const {
		index_type i = iSell / (nY+2);
		index_type j = iSell % (nY+2);
		return idxEdgeCell(i,j,dir);
	}
	inline index_type idxNodeCell( index_type i, index_type j, index_type dir ) const {
		if(i == 0 || i == nX+1 || j == 0 || j == nY+1)
			return -1;

		switch(dir)	{
		case 0 : return iNode(i-1,j-1); 
		case 1 : return iNode(i  ,j-1); 
		case 2 : return iNode(i  ,j  ); 
		case 3 : return iNode(i-1,j  ); 
		default : return -1;
		}		
	}
	inline index_type idxNodeCell( index_type iSell, index_type dir ) const {
		index_type i = iSell / (nY+2);
		index_type j = iSell % (nY+2);
		return idxNodeCell(i,j,dir);		
	}
	inline void getFourNodesOfCell(index_type i, index_type j, index_type nodes[4])	{
		for(index_type l=0;l<4;l++)
			nodes[l] = idxNodeCell(i,j,l);
	}

	inline T getS(const Vector2<T> &_pos, const T *_phi)	{
		index_type r = (index_type)(_pos[0]*hInv+0.5f);
		index_type s = (index_type)(_pos[1]*hInv+0.5f);

		if ( r == nX+1 ) r--;
		if ( s == nY+1 ) s--;

		Vector2<T> corner = posS(r,s);
		T alpha = ( _pos[0] - corner[0] ) * hInv;
		T beta  = ( _pos[1] - corner[1] ) * hInv;

		return (1.f-alpha)*((1.f-beta)*_phi[iS(r  ,s)] + beta*_phi[iS(r  ,s+1)])
				   + alpha*((1.f-beta)*_phi[iS(r+1,s)] + beta*_phi[iS(r+1,s+1)]);
	}
	inline T getU(const Vector2<T> &_pos, const T *_u)	{
		index_type r = (index_type)(_pos[0]*hInv);
		index_type s = (index_type)(_pos[1]*hInv+0.5f);

		if ( r == nX   ) r--;
		if ( s == nY+1 ) s--;

		Vector2<T> corner = posU(r,s);
		T alpha = ( _pos[0] - corner[0] ) * hInv;
		T beta  = ( _pos[1] - corner[1] ) * hInv;

		return (1.f-alpha)*((1.f-beta)*_u[iU(r  ,s)] + beta*_u[iU(r  ,s+1)])
				   + alpha*((1.f-beta)*_u[iU(r+1,s)] + beta*_u[iU(r+1,s+1)]);
	}
	inline T getV(const Vector2<T> &_pos, const T *_v)	{
		index_type r = (index_type)(_pos[0]*hInv+0.5f);
		index_type s = (index_type)(_pos[1]*hInv);

		if ( r == nX+1 ) r--;
		if ( s == nY   ) s--;

		Vector2<T> corner = posV(r,s);
		T alpha = ( _pos[0] - corner[0] ) * hInv;
		T beta  = ( _pos[1] - corner[1] ) * hInv;

		return (1.f-alpha)*((1.f-beta)*_v[iV(r  ,s)] + beta*_v[iV(r  ,s+1)])
				   + alpha*((1.f-beta)*_v[iV(r+1,s)] + beta*_v[iV(r+1,s+1)]);
	}
	inline T getNode(const Vector2<T> &_pos, const T *_phi)	{
		index_type r = (index_type)(_pos[0]*hInv);
		index_type s = (index_type)(_pos[1]*hInv);

		if ( r == nX ) r--;
		if ( s == nY ) s--;

		Vector2<T> corner = posNode(r,s);
		T alpha = ( _pos[0] - corner[0] ) * hInv;
		T beta  = ( _pos[1] - corner[1] ) * hInv;

		return (1.f-alpha)*((1.f-beta)*_phi[iNode(r  ,s)] + beta*_phi[iNode(r  ,s+1)])
				   + alpha*((1.f-beta)*_phi[iNode(r+1,s)] + beta*_phi[iNode(r+1,s+1)]);
	}
};

};		// geometryspace

#include "MacGrid2D_Functions.h"

#endif