#ifndef _FLUIDSPACE_CFD_INTERFACE_TRACKING_SRL2D_H_
#define _FLUIDSPACE_CFD_INTERFACE_TRACKING_SRL2D_H_

#include "../../../Principle/Geometry/SRL/SRL2D.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../FluidLibrary/InternalBoundary.h"
#include "InterfaceTracking.h"

namespace fluidspace {

using geometryspace::MacGrid2D;
using geometryspace::SRLNode2D;

class NATUREDLL InterfaceTracking_SRL2D : public InterfaceTracking {
public : 
	typedef geometryspace::SRL2D<real>::index_type index_type;

public : 
	InterfaceTracking_SRL2D() : srl(NULL), srl0(NULL) {}
	InterfaceTracking_SRL2D(index_type nX, index_type nY, real h, real narrowband = 4, bool uniform_quadrature = false) : srl(NULL), srl0(NULL) 
	{ set(nX, nY, h, narrowband, uniform_quadrature); }

	~InterfaceTracking_SRL2D() { release(); }

	void release();

	void set(index_type nX, index_type nY, real h, real narrowband = 5, bool uniform_quadrature = false);
	void setPhi(const geometryspace::DistanceField2D<real>& field) { srl0->setPhi(field); srl->setPhi(field); } 

	void addSource(const geometryspace::DistanceField2D<real>& field) { srl->addSource(field); }

	const index_type get_nX() const { return srl->phi.nX; }
	const index_type get_nY() const { return srl->phi.nY; }

	const real get_Lx() const { return srl->Lx; }
	const real get_Ly() const { return srl->Ly; }
	const real get_h() const { return srl->h; }

	const Vector2<real> get_offset2() const { return srl->phi.offset; } 

	const real getTotalVolume() const { if(isbdry.is_empty()) return srl->getTotalVolume(); else return srl->getTotalVolume(isbdry); }	
	const real getTotalVolume(const Grid2<bool,real>& isbdry) const { return srl->getTotalVolume(isbdry); }

	void getPhi(Grid2<real,real>&) const;	
	void getCurvature(Grid2<real,real>&) const;

	void update_InternalBoundary(const InternalBoundary_Implicit2D*);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void calc_narrow_region(real narrowband);
	template<typename D, int N> void add_narrowband_subcell_interface(Grid2<D,real>&, std::vector<SRLNode2D<real,N>>&, const D tag_narrow);

	template<typename D> void subtract_narrowband_kinematicbdry(Grid2<D,real>&, const D tag_far);

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void createDrawableInterface();
	void createDrawableInterface(const Grid2b& isbdry);

	void glDrawInterface(const virtualspace::Material& mat, bool filled) const;

	void glDrawPhi(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5, GLfloat ptsize = 3.0f) const;	
	void glDrawPhi_SubGrid(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5, GLfloat ptsize = 3.0f) const;	
	void glDrawGrid(const RGBA&, GLfloat width) const;
	void glDrawBoundingBox(const RGBA&, GLfloat width) const;
	void glDrawAuxiliary() const;

	const unsigned int getDimension() const { return 2;}
	const bool isEulerian() const { return true; }

	const bool load(const string& filename);
	const bool save(const string& filename) const;

	const bool load_compressed(const string& filename);
	const bool save_compressed(const string& filename) const;

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, MAX(srl->Lx,srl->Ly) * 1.42f); }
	const Vector3f viewRef() const { return Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, 0.0f);		}
	const Vector3f viewUp() const  { return Vector3f(0.0f, 1.0f, 0.0f);						}
	const float viewTheta() const  { return 45.0f;											}

	const std::string name_class() const { return XmlElementName(); }	
	static const std::string XmlElementName() { return std::string("IntefaceTracking_SRL2D"); }

public : 
	Grid2<bool,real>			isbdry;
	geometryspace::SRL2D<real> *srl, *srl0;

	// variable for member function _ InterfaceTracking_SRL2D::add_narrowband_subcell_interface()
	std::vector<int>			narrow_region;		
};

// Need to Parallelize SRL
template<typename D, int N> 
void InterfaceTracking_SRL2D::add_narrowband_subcell_interface(Grid2<D,real>& tag, std::vector<SRLNode2D<real,N>>& nodes, const D tag_narrow) {
	
#ifdef USE_OPENMP	
	int my_end = int(nodes.size());
	int i;

#pragma omp parallel for
#else	
	std::size_t my_end = nodes.size();
	std::size_t i;
#endif
	for( i = 0; i < my_end; ++i ) {			
		if(nodes[i].has_coarse_interface() == false && nodes[i].has_interface() == true) {

			index_type ii = nodes[i].index[0], jj = nodes[i].index[1];

			for(std::size_t k=0;k<narrow_region.size();k+=2) {
				index_type p = ii + narrow_region[k  ];
				index_type q = jj + narrow_region[k+1];
				
				if( 0 <= p && p < tag.nX-1 
				 && 0 <= q && q < tag.nY-1 )
					tag(p, q) = tag_narrow;
			}
		}
	}
}

template<typename D> 
void InterfaceTracking_SRL2D::subtract_narrowband_kinematicbdry(Grid2<D,real>& tag, const D tag_far) {

	if(isbdry.is_empty())
		return;

	for(index_type i=0;i<tag.nX-1;++i) {
		for(index_type j=0;j<tag.nY-1;++j) {

			if(isbdry(i,j) && isbdry(i+1,j) && isbdry(i,j+1) && isbdry(i+1,j+1))
				tag(i,j) = tag_far;

		}
	}	
}

};		// fluidspace

#endif