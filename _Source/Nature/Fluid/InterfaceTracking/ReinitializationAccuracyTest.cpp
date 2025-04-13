#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/ReinitializationAccuracyTest.h>
#include <Principle/Geometry/Levelset/TVD_RungeKutta.h>
#include <Principle/Geometry/Levelset/Levelset_AccuracyTest.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

ReinitializationAccuracyTest::ReinitializationAccuracyTest() {	
	tracking = NULL;
	redist = NULL;
	redist_srl = NULL;
	ioHandler = new ReinitializationAccuracyTest_IOHandler(this);

	//for(size_t iter=0;iter<10;++iter) {
	//	real v[5], h[4], dx_m, dx_p;
	//	for(int i=0;i<5;++i) { v[i] = rand() * 0.0001f; cout << v[i] << " "; }
	//	for(int i=0;i<4;++i) h[i] = 1;
	//	cout << endl;

	//	eno2nd_adative(dx_m, dx_p, v, h);	cout << dx_m << " " << dx_p << endl;
	//	dx_m = eno2nd(v, h[0], true);
	//	dx_p = eno2nd(v, h[0], false);		cout << dx_m << " " << dx_p << endl;
	//}	
	//exit(1);

	//double * weight = new double[3];
	//for(int i=0;i<3;++i) {
	//	weight_1st_derivative(weight, 3-i, 3, 1.0);
	//	cout << weight[0] << " " << weight[1] << " " << weight[2] << endl;
	//}

	//double weno[3];	cout << endl;
	//weno_weight_from_eno(weno, 2, 1.0);	cout << weno[0] << " " << weno[1] << " " << weno[2] << endl;
	//weno_weight_from_eno(weno, 3, 1.0);	cout << weno[0] << " " << weno[1] << " " << weno[2] << endl;
}

void ReinitializationAccuracyTest::release() {
	if(tracking) delete tracking; tracking = NULL;
	if(redist) delete redist; redist = NULL;
	if(redist_srl) delete redist_srl; redist_srl = NULL;
}

void ReinitializationAccuracyTest::setExperiment_Russo() {

	bool print_max_error = false;
	real band = 5; real cfl = 0.5;
	size_t min_iter = 0;

	size_t res_base = 64;
	size_t icr_res_step = 1;

	release();

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field
	real L = 4;
	RussoDistanceField2D<real> russo(Vector2<real>(L*0.5, L*0.5));	
	CircleDistanceField2D<real> circle(Vector2<real>(L*0.5, L*0.5), 1.0);	

	//////////////////////////////////////////////////////////////////////////
	// Setting Drawing Options
	dBoundingBox = dLiquid = true;
	dGrid = dPhi = dAuxiliary = false;	

	//////////////////////////////////////////////////////////////////////////
	// Structure
	InterfaceTracking_Levelset2D* lv2D = new InterfaceTracking_Levelset2D(10, 10, 1);
	tracking = lv2D;

	InterfaceTracking_SRL2D* srl2D = new InterfaceTracking_SRL2D();
	tracking = srl2D;

	Levelset_InfNorm<real> infnorm(2.0);
	Levelset_L1Norm<real> l1norm(2.0);

	real error_lv[4][10], error_srl[4][10];
	real max_error_phi_lv[4][10], max_error_phi_srl[4][10];
	real max_error_exact_lv[4][10], max_error_exact_srl[4][10];

	for(size_t i=0;i<icr_res_step;++i) {
		grid_index_type res = res_base * pow(float(2),float(i));
		real h = L / real(res);

		//redist = new Levelset_RedistancingFMM<real>(band);
		//redist = new Levelset_Redistancing_PDE_1stUpwind<real, RungeKutta_Euler>(band, cfl, min_iter);
		//redist = new Levelset_Redistancing_PDE_2ndUpwind<real, RungeKutta_RK2>(band, cfl, min_iter);
		redist = new Levelset_Redistancing_PDE_ENO<real, RungeKutta_RK3>(band, cfl, min_iter);
		//redist = new Levelset_Redistancing_PDE_WENO<real, RungeKutta_RK3>(band, cfl, min_iter);

		//redist_srl = new SRL_Redistancing_PDE_1stUpwind(band, 0.49f, true, min_iter);
		//redist_srl = new SRL_Redistancing_PDE_2ndUpwind(band, cfl, true, min_iter);
		redist_srl = new SRL_Redistancing_PDE_WENO3rd(band, cfl, true, min_iter);
		//redist_srl = new SRL_Redistancing_PDE_WENO(band, cfl, true, min_iter);
		//redist_srl = new SRL_Redistancing_Direct(band);

		lv2D->set(res, res, h);
		lv2D->setPhi(russo);	
		//lv2D->setPhi(circle);	

		srl2D->set(res, res, h, band);
		srl2D->setPhi(russo);
		//srl2D->setPhi(circle);	

		cout << "Redistancing " << res << " Grid" << endl;
		if(redist) {
			redist->redistance(*lv2D);
		}
		if(redist_srl) {
			redist_srl->redistance(*srl2D);
		}
		cout << "Redistancing " << res << " Grid Done" << endl << endl;

		real narrowband = 0; 

		narrowband = 1.4;	infnorm.set_narrowband(narrowband); l1norm.set_narrowband(narrowband);
		error_lv[0][i] = infnorm.norm(*lv2D, circle); max_error_phi_lv[0][i] = infnorm.max_error_phi; max_error_exact_lv[0][i] = infnorm.max_error_exact;
		error_lv[1][i] = l1norm.norm(*lv2D, circle);		

		narrowband = 3.0;	infnorm.set_narrowband(narrowband); l1norm.set_narrowband(narrowband);
		error_lv[2][i] = infnorm.norm(*lv2D, circle); max_error_phi_lv[2][i] = infnorm.max_error_phi; max_error_exact_lv[2][i] = infnorm.max_error_exact;
		error_lv[3][i] = l1norm.norm(*lv2D, circle);		

		narrowband = 1.4;	infnorm.set_narrowband(narrowband); l1norm.set_narrowband(narrowband);
		error_srl[0][i] = infnorm.norm(*srl2D->srl, circle); max_error_phi_srl[0][i] = infnorm.max_error_phi; max_error_exact_srl[0][i] = infnorm.max_error_exact;
		error_srl[1][i] = l1norm.norm(*srl2D->srl, circle);

		narrowband = 3.0;	infnorm.set_narrowband(narrowband); l1norm.set_narrowband(narrowband);
		error_srl[2][i] = infnorm.norm(*srl2D->srl, circle); max_error_phi_srl[2][i] = infnorm.max_error_phi; max_error_exact_srl[2][i] = infnorm.max_error_exact;
		error_srl[3][i] = l1norm.norm(*srl2D->srl, circle);

		if(redist) delete redist; redist = NULL;
		if(redist_srl) delete redist_srl; redist_srl = NULL;
	}

	for(size_t i=0;i<4;++i) {
		cout << i << "th error_lv Analysis" << endl;
		cout << error_lv[i][0] << endl;

		for(size_t step=1;step<icr_res_step;++step) {
			cout << error_lv[i][step] << " " << error_lv[i][step-1] / error_lv[i][step] << endl;
		}
		cout << endl;
	}
	
	if(print_max_error) {
		for(size_t i=0;i<4;++i) {		
			cout << i << "th error_lv Max Error" << endl;
			for(size_t step=0;step<icr_res_step;++step) {
				if(i % 2 == 0)
					cout << max_error_exact_lv[i][step] << " " << max_error_phi_lv[i][step] << endl;
			}

			cout << endl;
		}
	}	

	for(size_t i=0;i<4;++i) {
		cout << i << "th error_srl Analysis" << endl;
		cout << error_srl[i][0] << endl;
		for(size_t step=1;step<icr_res_step;++step) {
			cout << error_srl[i][step] << " " << error_srl[i][step-1] / error_srl[i][step] << endl;
		}
		cout << endl;
	}

	if(print_max_error) {
		for(size_t i=0;i<4;++i) {		
			cout << i << "th error_srl Max Error" << endl;
			for(size_t step=0;step<icr_res_step;++step) {
				if(i % 2 == 0)
					cout << max_error_exact_srl[i][step] << " " << max_error_phi_srl[i][step] << endl;
			}
			cout << endl;
		}
	}

	tracking->createDrawableInterface();	
}

void ReinitializationAccuracyTest::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	glDraw();	
}

void ReinitializationAccuracyTest::glDraw() const {

	if(tracking == NULL) return;

	Material mat;
	mat.setAmbient(RGBA::blue());
	mat.setDiffuse(RGBA::blue());
	mat.setSpecular(RGBA::blue() * 0.25f);
	mat.setShininess(8.5f);

	if(dLiquid)															{ tracking->glDrawInterface(mat, true);	} 
	if(dPhi && tracking->getDimension() == 2 && tracking->isEulerian())	{ tracking->glDrawPhi();				}

	if(dAuxiliary) {
		InterfaceTracking_SRL2D* srl = dynamic_cast<InterfaceTracking_SRL2D*>(tracking);
		if(srl != NULL && dPhi == false)	srl->glDrawAuxiliary();
		if(srl != NULL && dPhi == true)		srl->glDrawPhi_SubGrid();
	}

	if(dGrid && tracking->getDimension() == 2 && tracking->isEulerian())	{ 
		tracking->glDrawGrid();			
	}
	if(dBoundingBox && tracking->isEulerian())								{ 
		tracking->glDrawBoundingBox();	
	}
}

ReinitializationAccuracyTest_IOHandler::ReinitializationAccuracyTest_IOHandler(ReinitializationAccuracyTest *_s) {
	sim = _s;
}

ReinitializationAccuracyTest_IOHandler::~ReinitializationAccuracyTest_IOHandler() {
}

void ReinitializationAccuracyTest_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'b' : sim->dBoundingBox = !sim->dBoundingBox; break;
		case 'g' : sim->dGrid = !sim->dGrid; break;
		case 'l' : sim->dLiquid = !sim->dLiquid; break;		
		case 'p' : sim->dPhi = !sim->dPhi; break;
		case 'a' : sim->dAuxiliary = !sim->dAuxiliary; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void ReinitializationAccuracyTest_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void ReinitializationAccuracyTest_IOHandler::mouseMotionHandler(int x, int y) {
}
