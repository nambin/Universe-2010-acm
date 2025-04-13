#ifndef _FLUIDSPACE_CFD_INTERFACE_TRACKING_SLC3D_H_
#define _FLUIDSPACE_CFD_INTERFACE_TRACKING_SLC3D_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../SLC/slcSurface.H"
#include "InterfaceTracking.h"

namespace fluidspace {

using geometryspace::MacGrid3D;

class NATUREDLL InterfaceTracking_SLC3D : public InterfaceTracking, public geometryspace::MacGrid3D<double> {
public : 
	InterfaceTracking_SLC3D() : s(NULL), s0(NULL) {}
	InterfaceTracking_SLC3D(double lc, double uc, int max_lv) : s(NULL), s0(NULL) { set(lc, uc, max_lv); }

	virtual void set(double lc, double uc, int max_lv);
	virtual void setPhi(const geometryspace::DistanceField3D<real> &);

	const real get_Lx() const { return uc - lc; }
	const real get_Ly() const { return uc - lc; }
	const real get_Lz() const { return uc - lc; }
	const real get_h() const { return dx; }

	const real getTotalVolume() const { return -1; }

	void getPhi(Grid3<double,double>& _phi) const;

	void exportOBJ(const std::string& filename) const { char str[256]; strcpy_s(str, 256, filename.c_str()); s->objDump(str); }

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void createDrawableInterface() {}

	void glDrawInterface(const virtualspace::Material& mat, bool filled) const;
	void glDrawPhi(const RGBA&, const RGBA&, real band, GLfloat pointsize) const {}
	void glDrawGrid(const RGBA&, GLfloat linewidth) const {}
	void glDrawBoundingBox(const RGBA&, GLfloat width) const;

	const unsigned int getDimension() const { return 3; }
	const bool isEulerian() const { return true; }	

	const bool load(const string& filename);
	const bool save(const string& filename) const;

	const bool load_compressed(const string& filename)		 { return load(filename); }
	const bool save_compressed(const string& filename) const { return save(filename); }

	const std::string name_class() const { return XmlElementName(); }
	static const std::string XmlElementName() { return std::string("InterfaceTracking_SLC3D"); }

public : 
	double	lc, uc;
	int		max_level;

	slcspace::SlcSurface*	s0;
	slcspace::SlcSurface*	s;
};

};		// fluidspace

#endif