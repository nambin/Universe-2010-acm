#pragma once

#include "../../../Principle/Geometry/SRL/SRL3D.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../FluidLibrary/InternalBoundary.h"
#include "InterfaceTracking.h"

namespace fluidspace {

using geometryspace::MacGrid3D;
using geometryspace::SRLNode3D;

class NATUREDLL InterfaceTracking_SRL3D : public InterfaceTracking {
public : 
	typedef geometryspace::SRL3D<real>::index_type index_type;

public : 
	InterfaceTracking_SRL3D(); 
	InterfaceTracking_SRL3D(index_type nX, index_type nY, index_type nZ, real h, real narrowband = 4, bool uniform_quadrature = false) : srl(NULL), srl0(NULL) 
	{ set(nX, nY, nZ, h, narrowband, uniform_quadrature); }

	~InterfaceTracking_SRL3D() { release(); }

	void release();

	void set(index_type nX, index_type nY, index_type nZ, real h, real narrowband = 5, bool uniform_quadrature = false);
	void setPhi(const geometryspace::DistanceField3D<real>& field) { srl0->setPhi(field); (*srl) = (*srl0); } 

	void addSource(const geometryspace::DistanceField3D<real>& field) { srl->addSource(field); }

	void exportOBJ(const std::string& filename, const Transform<real>& trans) const { srl->mesh.exportOBJ(filename, trans); }

	const index_type get_nX() const { return srl->phi.nX; }
	const index_type get_nY() const { return srl->phi.nY; }
	const index_type get_nZ() const { return srl->phi.nZ; }

	const real get_Lx() const { return srl->Lx; }
	const real get_Ly() const { return srl->Ly; }
	const real get_Lz() const { return srl->Lz; }
	const real get_h() const { return srl->h; }

	const Vector3<real> get_offset3() const { return srl->phi.offset; } 

	const real getTotalVolume() const { if(isbdry.is_empty()) return srl->getTotalVolume(); else return srl->getTotalVolume(isbdry); }	
	const real getTotalVolume(const Grid3<bool,real>& isbdry) const { return srl->getTotalVolume(isbdry); }

	void getPhi(Grid3<real,real>& _phi) const;

	void getCurvature(Grid3<real,real>&) const;

	void update_InternalBoundary(const InternalBoundary_Implicit3D*);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void calc_narrow_region(real narrowband);
	template<typename D, int N> void add_narrowband_subcell_interface(Grid3<D,real>&, std::vector<SRLNode3D<real,N>>&, const D tag_narrow);

	template<typename D> void subtract_narrowband_kinematicbdry(Grid3<D,real>&, const D tag_far);

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void createDrawableInterface();
	void createDrawableInterface(const Grid3b&);

	void glDrawInterface(const virtualspace::Material& mat, bool filled) const;
	void glDrawInnerBoundary(const virtualspace::Material& mat) const;

	void glDrawPhi(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5, GLfloat ptsize = 3.0f) const;	
	void glDrawPhi_SubGrid(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5, GLfloat ptsize = 3.0f) const;	
	void glDrawGrid(const RGBA&, GLfloat width) const;
	void glDrawBoundingBox(const RGBA&, GLfloat width) const;
	void glDrawAuxiliary() const;

	const unsigned int getDimension() const { return 3;}
	const bool isEulerian() const { return true; }

	const bool load(const string& filename);
	const bool save(const string& filename) const;

	const bool load_compressed(const string& filename);
	const bool save_compressed(const string& filename) const;

	const std::string name_class() const { return XmlElementName(); }
	static const std::string XmlElementName() { return std::string("InterfaceTracking_SRL3D"); }

public : 
	Grid3<bool,real>			isbdry;
	geometryspace::SRL3D<real> *srl, *srl0;

	// variable for member function _ InterfaceTracking_SRL3D::add_narrowband_subcell_interface()
	std::vector<int>			narrow_region;	
};

// Need to Parallelize SRL
template<typename D, int N> 
void InterfaceTracking_SRL3D::add_narrowband_subcell_interface(Grid3<D,real>& tag, std::vector<SRLNode3D<real,N>>& nodes, const D tag_narrow) {

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

			index_type ii = nodes[i].index[0], jj = nodes[i].index[1], kk = nodes[i].index[2];

			for(std::size_t k=0;k<narrow_region.size();k+=3) {
				index_type p = ii + narrow_region[k  ];
				index_type q = jj + narrow_region[k+1];
				index_type r = kk + narrow_region[k+2];

				if( 0 <= p && p < tag.nX-1 
				 && 0 <= q && q < tag.nY-1 
				 && 0 <= r && r < tag.nZ-1 )
					tag(p,q,r) = tag_narrow;

			}
		}
	}
}

template<typename D> 
void InterfaceTracking_SRL3D::subtract_narrowband_kinematicbdry(Grid3<D,real>& tag, const D tag_far) {

	if(isbdry.is_empty())
		return;

	for(index_type i=0;i<tag.nX-1;++i) {
		for(index_type j=0;j<tag.nY-1;++j) {
			for(index_type k=0;k<tag.nZ-1;++k) {

				if(isbdry(i,j,k  ) && isbdry(i+1,j,k  ) && isbdry(i,j+1,k  ) && isbdry(i+1,j+1,k  )
				&& isbdry(i,j,k+1) && isbdry(i+1,j,k+1) && isbdry(i,j+1,k+1) && isbdry(i+1,j+1,k+1)	)

					tag(i,j,k) = tag_far;

			}
		}
	}	
}

};			// fluidspace