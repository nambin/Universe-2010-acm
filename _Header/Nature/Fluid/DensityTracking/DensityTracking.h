#ifndef _FLUIDSPACE_CFD_DENSITY_TRACKING_H_
#define _FLUIDSPACE_CFD_DENSITY_TRACKING_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../../../Principle/Geometry/MacGrid/MacGrid2D.h"
#include "../../../Principle/Geometry/MacGrid/MacGrid3D.h"
#include "../../../Principle/Geometry/MacGrid/MacGridLibrary.h"
#include "../../../Principle/Geometry/DistanceFieldLibrary/DistanceFields.h"
#include "../../../Illumination/RealtimeShader/cgShader/cgShader_VolumeRendering.h"

namespace fluidspace {

using geometryspace::MacGrid2D;
using geometryspace::MacGrid3D;
using geometryspace::DistanceField2D;
using geometryspace::DistanceField3D;
using shaderspace::cgShader_VolumeRendering;

class NATUREDLL DensityTracking {
public : 
	DensityTracking() {}

	virtual void post_process() const { log_info << "DensityTracking _ Have no Post Process" << std::endl; }

	virtual const real getTotalDensity() const = 0;

	virtual const unsigned int getDimension() const = 0;
	virtual void glDraw(const RGBA& c = RGBA::white()) const = 0;
};

class NATUREDLL DensityTracking2D : public DensityTracking, public MacGrid2D<real> {
public : 
	DensityTracking2D() {}

	void set(grid_index_type nX, grid_index_type nY, real h);
	void set(const Grid2r&);
	void set(const DistanceField2D<real>&);

	const real getTotalDensity() const;

	const unsigned int getDimension() const { return 2; }

	void glDraw(const RGBA& c = RGBA::white()) const;
	static void glDraw(const Grid2r&, const RGBA& c = RGBA::white());

public : 
	Grid2r density, density0;
};

class NATUREDLL DensityTracking3D : public DensityTracking, public MacGrid3D<real> {
public : 
	DensityTracking3D() { shader = NULL; isUpdated = false; }
	~DensityTracking3D() { if(shader) delete shader; }

	//////////////////////////////////////////////////////////////////////
	// Setter
	void set(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h);
	void set(const Grid3r&);
	void set(const DistanceField3D<real>&);

	const real getTotalDensity() const;

	const unsigned int getDimension() const { return 3; }

	//////////////////////////////////////////////////////////////////////
	// Getter
	const Grid3r& getDensity() const  { return density;	 }
	const Grid3r& getDensity0()	const { return density0; }
	Grid3r&	getDensity()	{ isUpdated = true; return density;	}
	Grid3r& getDensity0()	{ isUpdated = true; return density0;}

	//////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw(const RGBA& c = RGBA::white()) const;	
	void enableDrawBoundingBox() const  { if(shader) shader->dBoundingBox = true;  	}
	void disableDrawBoundingBox() const { if(shader) shader->dBoundingBox = false;	}

	void initialize_Shader(unsigned int width, unsigned int height);	

protected : 
	Grid3r			density, density0;
	mutable bool	isUpdated;

	// For Rendering
	cgShader_VolumeRendering*	shader;
};

};		// fluidspace

#endif