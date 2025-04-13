#ifndef _FLUIDSPACE_CFD_INTERFACE_TRACKING_LEVELSET3D_H_
#define _FLUIDSPACE_CFD_INTERFACE_TRACKING_LEVELSET3D_H_

#include "../../../Principle/Geometry/Levelset/Levelset3D.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "InterfaceTracking.h"

namespace fluidspace {

using geometryspace::MacGrid3D;

class NATUREDLL InterfaceTracking_Levelset3D : public InterfaceTracking, public geometryspace::Levelset3D<real> {
public : 
	typedef geometryspace::Levelset3D<real>::index_type index_type;

public : 
	InterfaceTracking_Levelset3D() {}
	InterfaceTracking_Levelset3D(index_type nX, index_type nY, index_type nZ, real h) { set(nX, nY, nZ, h); }

	virtual void set(index_type nX, index_type nY, index_type nZ, real h);
	virtual void setPhi(const Grid3r&);
	virtual void setPhi(const geometryspace::DistanceField3D<real> &);

	const index_type get_nX() const { return phi.nX; }
	const index_type get_nY() const { return phi.nY; }
	const index_type get_nZ() const { return phi.nZ; }

	const real get_Lx() const { return Lx; }
	const real get_Ly() const { return Ly; }
	const real get_Lz() const { return Lz; }
	const real get_h() const { return h; }

	const Vector3<real> get_offset3() const { return phi.offset; } 

	const real getTotalVolume() const { return geometryspace::Levelset3D<real>::getTotalVolume(); }
	const real getTotalVolume(const Grid3<bool,real>& isbdry) const { return geometryspace::Levelset3D<real>::getTotalVolume(isbdry); }
	
	void getPhi(Grid3<real,real>& _phi) const;
	void getNormal(Grid3<Vector3<real>,real>& _normal) const;
	void getCurvature(Grid3<real,real>& _curvature) const;

	const bool load_compressed(const string& filename);
	const bool save_compressed(const string& filename) const;

	void exportOBJ(const std::string& filename, const Transform<real>& trans) const { mesh.exportOBJ(filename, trans); }

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void createDrawableInterface();
	void createDrawableInterface(const Grid3b& isbdry);
	
	void glDrawInterface(const virtualspace::Material& mat, bool filled) const;
	void glDrawPhi(const RGBA&, const RGBA&, real band, GLfloat pointsize) const {}
	void glDrawGrid(const RGBA&, GLfloat linewidth) const {}
	void glDrawBoundingBox(const RGBA&, GLfloat width) const;

	const unsigned int getDimension() const { return 3; }
	const bool isEulerian() const { return true; }	

	const std::string name_class() const { return XmlElementName(); }
	static const std::string XmlElementName() { return std::string("InterfaceTracking_Levelset3D"); }

public : 
	Grid3r	phi0;	
};

};		// fluidspace

#endif