#ifndef _FLUIDSPACE_INTERFACE_TRACKING_TEST_H_
#define _FLUIDSPACE_INTERFACE_TRACKING_TEST_H_

#include "../CFDsolver/CFDsolver2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL3D.h"
#include "../InterfaceTracking/InterfaceTracking_SLC3D.h"
#include "../InterfaceTracking/InterfaceIntegrator_Levelset.h"
#include "../InterfaceTracking/InterfaceIntegrator_SRL.h"
#include "../InterfaceTracking/InterfaceIntegrator_SLC.h"
#include "../FluidScene/FluidScene.h"
#include "../CFDsolver/CFDIntegrator.h"
#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"

namespace fluidspace {

using arithmetic::Vector3f;

class InterfaceTrackingTest;
class NATUREDLL InterfaceTrackingTest_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	InterfaceTrackingTest_IOHandler(InterfaceTrackingTest *_sim);
	~InterfaceTrackingTest_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	InterfaceTrackingTest *sim;
};

class NATUREDLL InterfaceTrackingTest : public FluidSimulator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	InterfaceTrackingTest();
	~InterfaceTrackingTest() { release(); }

	void post_process() const;

	void release();
	
	//////////////////////////////////////////////////////////////////////////
	// Initialization	
	template<typename FIELD> void set_InterfaceTracking(grid_index_type nX, grid_index_type nY, real h, const FIELD&);
	template<typename FIELD> void set_InterfaceTracking(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, const FIELD&);
	void set_InterfaceIntegrator();

	void setExperiment_Translation_2D();	
	void setExperiment_Zalesak_2D();

	void setExperiment_Translation_3D();
	void setExperiment_Zalesak_3D();
	void setExperiment_Enright_3D();
	void setExperiment_Object_3D(const std::string& obj_file);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	void saveBinary(std::size_t frame) const;

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt);

	const unsigned int getDimension() const { return -1; }

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { if(tracking) return tracking->viewEye();	return Vector3f(0,0,1); }
	const Vector3f viewRef() const { if(tracking) return tracking->viewRef();	return Vector3f(0,0,0); }
	const Vector3f viewUp() const  { if(tracking) return tracking->viewUp();	return Vector3f(0,1,0); }
	const float viewTheta() const  { if(tracking) return tracking->viewTheta();	return 45.0f;			}

public :
	VelocityField<real>*	vel;

	InterfaceTracking*		tracking;
	InterfaceIntegrator*	integrator;

	std::string				test_name;

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dLiquid;
	bool dGrid;	
	bool dPhi;
	bool dAuxiliary;
	bool dFilledMesh;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend InterfaceTrackingTest_IOHandler;
};

const real narrowband_srl = 4;
const bool uniform_quadrature = false;

const bool use_srl = true;
const bool use_slc = false;
const bool use_pls = false;

template<typename FIELD> 
void InterfaceTrackingTest::set_InterfaceTracking(grid_index_type nX, grid_index_type nY, real h, const FIELD& field) {

	if(use_srl) {
		InterfaceTracking_SRL2D* srl2D = new InterfaceTracking_SRL2D();
		srl2D->set(nX, nY, h, narrowband_srl, uniform_quadrature);
		srl2D->setPhi(field);
		tracking = srl2D;
	}	
	else {
		InterfaceTracking_Levelset2D* lv2D = new InterfaceTracking_Levelset2D(nX, nY, h);
		lv2D->setPhi(field);
		tracking = lv2D;
	}	

	tracking->createDrawableInterface();	
}

template<typename FIELD> 
void InterfaceTrackingTest::set_InterfaceTracking(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, const FIELD& field) {

	if(use_srl) {
		InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
		srl3D->set(nX, nY, nZ, h, narrowband_srl, uniform_quadrature);
		srl3D->setPhi(field);
		tracking = srl3D;
	}
	else if(use_slc) {
		int max_lv = 0; grid_index_type tmp = nX;
		while(true) { if(tmp == 1) break; tmp /= 2; max_lv++; }

		InterfaceTracking_SLC3D* slc3D = new InterfaceTracking_SLC3D(0, nX*h, max_lv);
		slc3D->setPhi(field);
		tracking = slc3D;
	}
	else {
		InterfaceTracking_Levelset3D* lv3D = NULL;

		if(!use_pls)	lv3D = new InterfaceTracking_Levelset3D(nX, nY, nZ, h);
		else			lv3D = new InterfaceTracking_PLS3D(nX, nY, nZ, h);

		lv3D->setPhi(field);
		tracking = lv3D;
	}	

	tracking->createDrawableInterface();	
}

inline void InterfaceTrackingTest::set_InterfaceIntegrator() {

	if(use_srl) {
		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(true));
		//itgSRL->set_Redistancing(new SRL_Redistancing_PDE_1stUpwind(narrowband_srl, 0.5, true, 0));
		//itgSRL->set_Redistancing(new SRL_Redistancing_PDE_2ndUpwind(narrowband_srl, 0.3, true, 0));
		itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO3rd(narrowband_srl, 0.3, true, 0));
		//itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.3, true, 0));
		//itgSRL->set_Redistancing(new SRL_Redistancing_Direct(narrowband_srl));
		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
		integrator = itgSRL;
	}
	else if(use_slc) {
		InterfaceIntegrator_SLC* itgSLC = new InterfaceIntegrator_SLC();
		integrator = itgSLC;
	}
	else {
		real narrow_lv = 5;

		if(!use_pls) {
			InterfaceIntegrator_Levelset* itgLv = new InterfaceIntegrator_Levelset();

			//itgLv->set_Advection(new CFDAdvection_SemiLagrangian<CFDTracking_RK2>());		
			itgLv->set_Advection(new CFDAdvection_BFECC(false));
			//itgLv->set_Advection(new CFDAdvection_2ndUpwind<RungeKutta_RK3>());
			//itgLv->set_Advection(new CFDAdvection_WENO<RungeKutta_RK3>());

			//itgLv->set_Redistancing(new Levelset_RedistancingFMM<real>(narrow_lv));
			//itgLv->set_Redistancing(new Levelset_Redistancing_PDE_1stUpwind<real>(narrow_lv, 0.5));
			//itgLv->set_Redistancing(new Levelset_Redistancing_PDE_2ndUpwind<real>(narrow_lv, 0.5));
			//itgLv->set_Redistancing(new Levelset_Redistancing_PDE_WENO<real>(narrow_lv, 0.5));
			itgLv->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
			integrator = itgLv;
		}
		else {
			InterfaceIntegrator_PLS* itgLv = new InterfaceIntegrator_PLS();
			itgLv->set_Advection(new CFDAdvection_SemiLagrangian<CFDTracking_EE>());	
			itgLv->set_Redistancing(new Levelset_RedistancingFMM<real>(narrow_lv));
			itgLv->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
			integrator = itgLv;
		}		
	}

	if(integrator)
		integrator->initialize(*tracking);
}

};		// fluidspace

#endif