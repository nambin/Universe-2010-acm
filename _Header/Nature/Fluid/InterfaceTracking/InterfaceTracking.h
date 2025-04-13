#pragma once

#include "../../../Principle/VirtualSpace/Material/Material.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid2.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3.h"
#include "../../Nature.h"

namespace fluidspace {

class NATUREDLL InterfaceTracking {
public : 
	virtual const std::string name_class() const = 0;

	InterfaceTracking();
	virtual ~InterfaceTracking() {}

	virtual void post_process() const { log_info << "InterfaceTracking _ Have no Post Process" << std::endl; }

	virtual void createDrawableInterface() = 0;
	virtual void createDrawableInterface(const Grid2b& isbdry) { basis::throwError("Not Implemented _ InterfaceTracking::createDrawableInterface(const Grid2b&)"); }
	virtual void createDrawableInterface(const Grid3b& isbdry) { basis::throwError("Not Implemented _ InterfaceTracking::createDrawableInterface(const Grid3b&)"); }

	virtual void exportOBJ(const std::string& filename, const Transform<real>& trans = Transform<real>()) const { basis::throwError("Not Implemented _ InterfaceTracking::exportOBJ()"); }

	virtual void glDrawInterface(const virtualspace::Material& mat = RGBA::blue(), bool filled = true) const = 0;
	virtual void glDrawPhi(const RGBA& cp = RGBA::red(), const RGBA &cn = RGBA::white(), real band = 5.0f, GLfloat ptsize = 3.0f) const = 0;	
	virtual void glDrawGrid(const RGBA& c = RGBA::white() * 0.7f, GLfloat linewidth = 0.5f) const = 0;
	virtual void glDrawBoundingBox(const RGBA& c = RGBA::white(), GLfloat linewidth = 1.5f) const = 0;
	virtual void glDrawAuxiliary() const {}

	virtual void glDrawInnerBoundary(const virtualspace::Material& mat = RGBA::brown()) const { throwError("Not Implemented _ InterfaceTracking::glDrawInnerBoundary(const virtualspace::Material&)"); }

	virtual const unsigned int getDimension() const = 0;
	virtual const bool isEulerian() const = 0;

	virtual const grid_index_type get_nX() const { basis::throwError("Not Implemented _ InterfaceTracking::get_nX()"); return 0; }
	virtual const grid_index_type get_nY() const { basis::throwError("Not Implemented _ InterfaceTracking::get_nY()"); return 0; }
	virtual const grid_index_type get_nZ() const { basis::throwError("Not Implemented _ InterfaceTracking::get_nZ()"); return 0; }

	virtual const real get_Lx() const { basis::throwError("Not Implemented _ InterfaceTracking::get_Lx()"); return 0; }
	virtual const real get_Ly() const { basis::throwError("Not Implemented _ InterfaceTracking::get_Ly()"); return 0; }
	virtual const real get_Lz() const { basis::throwError("Not Implemented _ InterfaceTracking::get_Lz()"); return 0; }
	virtual const real get_h() const { basis::throwError("Not Implemented _ InterfaceTracking::get_h()"); return 0; }

	virtual const Vector2<real> get_offset2() const { basis::throwError("Not Implemented _ InterfaceTracking::get_offset2()"); return Vector2<real>(0); } 
	virtual const Vector3<real> get_offset3() const { basis::throwError("Not Implemented _ InterfaceTracking::get_offset3()"); return Vector3<real>(0); } 

	virtual const real getTotalVolume() const { basis::throwError("Not Implemented _ InterfaceTracking::getTotalVolume()"); return 0; }
	virtual const real getTotalVolume(const Grid2<bool,real>& isbdry) const { basis::throwError("Not Implemented _ InterfaceTracking::getTotalVolume()"); return 0; }
	virtual const real getTotalVolume(const Grid3<bool,real>& isbdry) const { basis::throwError("Not Implemented _ InterfaceTracking::getTotalVolume()"); return 0; }

	virtual void getNormal(arithmetic::Grid2<Vector2<real>,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getNormal()"); }
	virtual void getNormal(arithmetic::Grid3<Vector3<real>,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getNormal()"); }

	virtual void getCurvature(arithmetic::Grid2<real,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getCurvature()"); }
	virtual void getCurvature(arithmetic::Grid3<real,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getCurvature()"); }

	virtual void getPhi(arithmetic::Grid2<real,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getPhi()"); }
	virtual void getPhi(arithmetic::Grid3<real,real>&) const { basis::throwError("Not Implemented _ InterfaceTracking::getPhi()"); }

	virtual void setPhi(const geometryspace::DistanceField2D<real>&) { basis::throwError("Not Implemented _ InterfaceTracking::setPhi()"); }
	virtual void setPhi(const geometryspace::DistanceField3D<real>&) { basis::throwError("Not Implemented _ InterfaceTracking::setPhi()"); }

	virtual void addSource(const geometryspace::DistanceField2D<real>&) { basis::throwError("Not Implemented _ InterfaceTracking::addSource()"); }
	virtual void addSource(const geometryspace::DistanceField3D<real>&) { basis::throwError("Not Implemented _ InterfaceTracking::addSource()"); }

	virtual const bool load(const string& filename) { basis::throwError("Not Implemented _ InterfaceTracking::load()"); return false; }
	virtual const bool save(const string& filename) const { basis::throwError("Not Implemented _ InterfaceTracking::save()"); return false; }

	virtual const bool load_compressed(const string& filename) { basis::throwError("Not Implemented _ InterfaceTracking::load_compressed()"); return false; }
	virtual const bool save_compressed(const string& filename) const { basis::throwError("Not Implemented _ InterfaceTracking::save_compressed()"); return false; }

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	virtual const Vector3f viewEye() const { return eye;	}
	virtual const Vector3f viewRef() const { return ref;	}
	virtual const Vector3f viewUp() const  { return up;		}
	virtual const float viewTheta() const  { return theta;	}

public : 
	Vector3f eye, ref, up;
	float theta;
};	
	
};		// fluidspace