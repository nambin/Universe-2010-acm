#ifndef _FLUIDSPACE_CFD_INTERFACE_TRACKING_LEVELSET2D_H_
#define _FLUIDSPACE_CFD_INTERFACE_TRACKING_LEVELSET2D_H_

#include "../../../Principle/Geometry/Levelset/Levelset2D.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "InterfaceTracking.h"

namespace fluidspace {

using geometryspace::MacGrid2D;

class NATUREDLL InterfaceTracking_Levelset2D : public InterfaceTracking, public geometryspace::Levelset2D<real> {
public : 
	typedef geometryspace::Levelset2D<real>::index_type index_type;

public : 
	InterfaceTracking_Levelset2D() {}
	InterfaceTracking_Levelset2D(index_type nX, index_type nY, real h) { set(nX, nY, h); }

	virtual void set(index_type nX, index_type nY, real h);
	virtual void setPhi(const Grid2r&);
	virtual void setPhi(const geometryspace::DistanceField2D<real> &);

	const index_type get_nX() const { return phi.nX; }
	const index_type get_nY() const { return phi.nY; }

	const real get_Lx() const { return Lx; }
	const real get_Ly() const { return Ly; }
	const real get_h() const { return h; }

	const Vector2<real> get_offset2() const { return phi.offset; } 

	const real getTotalVolume() const { return geometryspace::Levelset2D<real>::getTotalVolume(); }
	const real getTotalVolume(const Grid2<bool,real>& isbdry) const { return geometryspace::Levelset2D<real>::getTotalVolume(isbdry); }

	void getPhi(Grid2<real,real>& _phi) const;
	void getNormal(Grid2<Vector2<real>,real>& _normal) const;
	void getCurvature(Grid2<real,real>& _curvature) const;

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void createDrawableInterface();
	void createDrawableInterface(const Grid2b& isbdry);

	void glDrawInterface(const virtualspace::Material& mat = RGBA::blue(), bool filled = true) const;
	void glDrawPhi(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5.0f, GLfloat ptsize = 3.0f) const;
	void glDrawGrid(const RGBA& c = RGBA::white() * 0.3f, GLfloat linewidth = 0.5f) const;
	void glDrawBoundingBox(const RGBA& c = RGBA::white(), GLfloat linewidth = 1.5f) const;

	const unsigned int getDimension() const { return 2; }
	const bool isEulerian() const { return true; }

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f*Lx, 0.5f*Ly, MAX(Lx,Ly) * 1.42f); }
	const Vector3f viewRef() const { return Vector3f(0.5f*Lx, 0.5f*Ly, 0.0f);				}
	const Vector3f viewUp() const  { return Vector3f(0.0f, 1.0f, 0.0f);						}
	const float viewTheta() const  { return 45.0f;											}

	const std::string name_class() const { return XmlElementName(); }
	static const std::string XmlElementName() { return std::string("InterfaceTracking_Levelset2D"); }

public : 
	Grid2r	phi0;
};

};		// fluidspace

#endif