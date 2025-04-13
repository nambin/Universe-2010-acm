#pragma once

#include "SRL_FileIO.h"

namespace geometryspace {

template<typename T, int N>
const T SRL2D<T,N>::getTotalVolume() const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			vol += MacGridLibrary<T>::heaviside(phi(i,j) * hInv, static_cast<T>(0.0), static_cast<T>(1.0));;
	return vol * h * h;
}

template<typename T, int N>
const T SRL2D<T,N>::getTotalVolume(const Grid2<bool,T>& isbdry) const {
	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			if(isbdry(i,j) == false)
				vol += MacGridLibrary<T>::heaviside(phi(i,j) * hInv, static_cast<T>(0.0), static_cast<T>(1.0));;
	return vol * h * h;
}

template<typename T, int N>
void SRL2D<T,N>::createRenderingPrimitive() {
	SRL<real>::SRL_Interpolator interp(true);
	SRL<real>::SRL_Lerp			lerp;
	SRL<real>::SRL_Coarse_Lerp	coarse_lerp;

	Grid2<T,T> phiNode(nX * (N-1) + 1, nY * (N-1) + 1, h / T(N-1), h / T(N-1), Vector2<T>(0));
	for(index_type i=0;i<phiNode.nX;++i) {
		for(index_type j=0;j<phiNode.nY;++j) {
			index_type ii, jj;	Vector2<T> x;
			get_idx_normalized_coord(ii, jj, x, phiNode.pos(i,j));

			if(is_exist_node(ii,jj))	{ phiNode(i,j) = lerp(*this, ii, jj, x);		}
			else						{ phiNode(i,j) = coarse_lerp(*this, ii, jj, x);	}
		}
	}

	Levelset2D<T>::marchingTriangle(renderingTriangles, phiNode, Grid2<bool,T>());
}

template<typename T, int N>
void SRL2D<T,N>::createRenderingPrimitive(const Grid2<bool,T>& isbdry) {
//	createRenderingPrimitive(); return;

	if(!isbdry.is_empty()) {
		marchingTriangle(renderingTriangles, isbdry);
		marchingInterface(renderingSegments, isbdry);
	}
	else {
		Grid2<bool,T> _isbdry; _isbdry.resize(phi);
		_isbdry.set(false);
		_isbdry.set_boundary(true);
		marchingTriangle(renderingTriangles, _isbdry);
		marchingInterface(renderingSegments, _isbdry);
	}	
}

template<typename T, int N>
void SRL2D<T,N>::marchingTriangle(vector<Triangle2D<T>>& triangles, const Grid2<bool,T>& isbdry) {
	triangles.clear();

	vector<Triangle2D<T>>	tmp;
	Vector2<T>				org, h;

	T _phi[4];

	for(index_type i=0;i<phi.nX-1;++i)	{
		for(index_type j=0;j<phi.nY-1;++j)	{

			if(isbdry(i,j) && isbdry(i+1,j) && isbdry(i,j+1) && isbdry(i+1,j+1))
				continue;

			if(indices(i,j) == NOT_REFINE) {
				if(!(phi(i,j) >  0 && phi(i+1,j) >  0 && phi(i,j+1) >  0 && phi(i+1,j+1) >  0) 
				&& !(phi(i,j) <= 0 && phi(i+1,j) <= 0 && phi(i,j+1) <= 0 && phi(i+1,j+1) <= 0) ) {
					throwError("Error in SRL2D<T,N>::marchingTriangle(vector<Triangle2D<T>>& triangles, const Grid2<bool,T>& isbdry)");
				}					

				if(phi(i,j) > 0)
					continue;

				_phi[0] = _phi[1] = _phi[2] = _phi[3] = -1;
				h[0] = phi.dx*0.5f; h[1] = phi.dy*0.5f;

				for(index_type p=0;p<2;++p) {
					for(index_type q=0;q<2;++q) {

						if(isbdry(i+p,j+q))
							continue;

						org = phi.pos(i,j) + Vector2<T>(h[0]*p, h[1]*q);						
						
						contouring(tmp, _phi, org, h);

						basis::vectorLibrary::append(triangles, tmp);

					}
				}				
				continue;
			}

			const int half_N = N/2;

			const SRLNode2D<T,N>& node = nodes[indices(i,j)];

			for(int p=0;p<N-1;++p) {
				for(int q=0;q<N-1;++q) {

					if		( p <  half_N && q <  half_N && isbdry(i  ,j  ) ) continue;
					else if	( p >= half_N && q <  half_N && isbdry(i+1,j  ) ) continue;
					else if	( p <  half_N && q >= half_N && isbdry(i  ,j+1) ) continue;
					else if	( p >= half_N && q >= half_N && isbdry(i+1,j+1) ) continue;

					org = pos(node, p, q);
					h[0] = h_sub[p]; h[1] = h_sub[q];

					_phi[0] = node(p  ,q  );	_phi[1] = node(p+1,q  );
					_phi[2] = node(p+1,q+1);	_phi[3] = node(p  ,q+1);

					contouring(tmp, _phi, org, h);

					basis::vectorLibrary::append(triangles, tmp);	

				}
			}
		}
	}	
}

template<typename T, int N>
void SRL2D<T,N>::marchingInterface(vector<Segment2D<T>>& segments, const Grid2<bool,T>& isbdry) {
	segments.clear();

	vector<Segment2D<T>>	tmp;
	Vector2<T>				org, h;

	T _phi[4];
	bool bdry_xm, bdry_xp, bdry_ym, bdry_yp;

	for(index_type i=0;i<phi.nX-1;++i)	{
		for(index_type j=0;j<phi.nY-1;++j)	{

			if(isbdry(i,j) && isbdry(i+1,j) && isbdry(i,j+1) && isbdry(i+1,j+1))
				continue;

			if(indices(i,j) == NOT_REFINE) {
				if(!(phi(i,j) >  0 && phi(i+1,j) >  0 && phi(i,j+1) >  0 && phi(i+1,j+1) >  0) 
				&& !(phi(i,j) <= 0 && phi(i+1,j) <= 0 && phi(i,j+1) <= 0 && phi(i+1,j+1) <= 0) ) {
					throwError("Error in SRL2D<T,N>::marchingTriangle(vector<Triangle2D<T>>& triangles, const Grid2<bool,T>& isbdry)");
				}					

				if(phi(i,j) > 0 || (!isbdry(i,j) && !isbdry(i+1,j) && !isbdry(i,j+1) && !isbdry(i+1,j+1)))
					continue;

				_phi[0] = _phi[1] = _phi[2] = _phi[3] = -1;
				h[0] = phi.dx*0.5f; h[1] = phi.dy*0.5f;

				for(index_type p=0;p<2;++p) {
					for(index_type q=0;q<2;++q) {

						if(isbdry(i+p,j+q))
							continue;

						bdry_xm = bdry_xp = bdry_ym = bdry_yp = false;

						if(p==0)	bdry_xp = isbdry(i+1,j+q);
						else 		bdry_xm = isbdry(i  ,j+q);

						if(q==0)	bdry_yp = isbdry(i+p,j+1);
						else		bdry_ym = isbdry(i+p,j  );

						org = phi.pos(i,j) + Vector2<T>(h[0]*p, h[1]*q);				
						contouring(tmp, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp);

						basis::vectorLibrary::append(segments, tmp);
					}
				}				
				continue;
			}

			const int half_N = N/2;

			const SRLNode2D<T,N>& node = nodes[indices(i,j)];

			for(int p=0;p<N-1;++p) {
				for(int q=0;q<N-1;++q) {

					if		( p <  half_N && q <  half_N && isbdry(i  ,j  ) ) continue;
					else if	( p >= half_N && q <  half_N && isbdry(i+1,j  ) ) continue;
					else if	( p <  half_N && q >= half_N && isbdry(i  ,j+1) ) continue;
					else if	( p >= half_N && q >= half_N && isbdry(i+1,j+1) ) continue;

					int i_ = (p < half_N) ? 0 : 1;
					int j_ = (q < half_N) ? 0 : 1;

					bdry_xm = bdry_xp = bdry_ym = bdry_yp = false;

					if(p==half_N-1)		bdry_xp = isbdry(i+1,j+j_);
					else if(p==half_N)	bdry_xm = isbdry(i  ,j+j_);

					if(q==half_N-1)		bdry_yp = isbdry(i+i_,j+1);
					else if(q==half_N)	bdry_ym = isbdry(i+i_,j  );

					org = pos(node, p, q);
					h[0] = h_sub[p]; h[1] = h_sub[q];

					_phi[0] = node(p  ,q  );	_phi[1] = node(p+1,q  );
					_phi[2] = node(p+1,q+1);	_phi[3] = node(p  ,q+1);

					contouring(tmp, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp);

					basis::vectorLibrary::append(segments, tmp);			
				}
			}
		}
	}	
}

template<typename T, int N>
void SRL2D<T,N>::glDrawLiquid(const RGBA& color) const {
	glColor4f(color[0], color[1], color[2], color[3]);

	glBegin(GL_TRIANGLES);
	for(size_t i=0;i<renderingTriangles.size();i++)	{
		glVertex3f(renderingTriangles[i].v[0][0], renderingTriangles[i].v[0][1], 0);
		glVertex3f(renderingTriangles[i].v[1][0], renderingTriangles[i].v[1][1], 0);
		glVertex3f(renderingTriangles[i].v[2][0], renderingTriangles[i].v[2][1], 0);
	}
	glEnd();
}

template<typename T, int N> 
void SRL2D<T,N>::glDrawInterfaceEdges(const RGBA &color, GLfloat width) const {
	glColor4f(color[0], color[1], color[2], color[3]);
	glLineWidth(width);

	glBegin(GL_LINES);
	for(std::size_t i=0;i<renderingSegments.size();++i) {
		glVertex3f(renderingSegments[i].v[0][0], renderingSegments[i].v[0][1], 0);
		glVertex3f(renderingSegments[i].v[1][0], renderingSegments[i].v[1][1], 0);
	}

	glEnd();
}

template<typename T, int N>
void SRL2D<T,N>::glDrawPhi(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const {

	glPointSize(ptsize);

	glBegin(GL_POINTS);
	for(index_type i=1;i<phi.nX-1;i++)	{
		for(index_type j=1;j<phi.nY-1;j++)	{

			T scale = (band - abs(phi(i,j)) * hInv) / band;

			if(scale < 0)			continue;
			else if(phi(i,j) > 0)	glColor4f(cp[0] * scale, cp[1] * scale, cp[2] * scale, cp[3]);
			else					glColor4f(cn[0] * scale, cn[1] * scale, cn[2] * scale, cn[3]);

			Vector2<T> pos = phi.pos(i,j);
			glVertex3f(pos[0], pos[1], 0.001f);
		}
	}	
	glEnd();
}

template<typename T, int N>
void SRL2D<T,N>::glDrawPhi_SubGrid(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const {

	glPointSize(ptsize);

	Vector2<T> pos0, p;
	for(size_t i=0;i<nodes.size();++i) {
		pos0 = phi.pos(nodes[i].index[0], nodes[i].index[1]);

		glBegin(GL_POINTS);
		for(index_type k=0;k<N;++k) {
			for(index_type l=0;l<N;++l) {
				T val = nodes[i].val[k][l];
				T scale = (band - abs(val) * hInv) / band;

				if(scale < 0)		continue;
				else if(val > 0)	glColor4f(cp[0] * scale, cp[1] * scale, cp[2] * scale, cp[3]);
				else				glColor4f(cn[0] * scale, cn[1] * scale, cn[2] * scale, cn[3]);
				
				p = pos0 + Vector2<T>(pos_sub[k], pos_sub[l]);
				glVertex3f(p[0], p[1], 0);
			}			
		}		
		glEnd();
	}
}

template<typename T, int N>
void SRL2D<T,N>::glDrawGrid(const RGBA& c, GLfloat width) const {
	phi.glDrawGrid(c, width);	
}

template<typename T, int N>
void SRL2D<T,N>::glDrawSubGrid(const RGBA& c, GLfloat ptsize) const {

	glPointSize(ptsize);
	glColor4f( c[0], c[1], c[2], c[3] );

	Vector2<T> pos0, p;
	for(size_t i=0;i<nodes.size();++i) {
		pos0 = phi.pos(nodes[i].index[0], nodes[i].index[1]);
		
		glBegin(GL_POINTS);
		for(index_type k=0;k<N;++k) {
			for(index_type l=0;l<N;++l) {
				p = pos0 + Vector2<T>(pos_sub[k], pos_sub[l]);
				glVertex3f(p[0], p[1], 0);
			}			
		}		
		glEnd();
	}
}

};		// geometryspace
