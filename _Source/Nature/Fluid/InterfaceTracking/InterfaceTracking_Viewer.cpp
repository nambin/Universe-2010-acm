#include <stdafx.h>
#include <iostream>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_Viewer.h>

using namespace std;
using namespace arithmetic;
using namespace fluidspace;

InterfaceTracking_Viewer::InterfaceTracking_Viewer() {
	tracking = NULL;
	filename = "";
	frame = 0;

	dBoundingBox = dLiquid = dInnerBdry = true;
	dMesh = dGrid = dPhi = dAuxiliary = false;	

	ioHandler = new InterfaceTracking_Viewer_IOHandler(this);
}

void InterfaceTracking_Viewer::exportOBJ(const string& filename) const {
	tracking->exportOBJ(filename);
	std::cout << filename.c_str() << " export obj is done" << std::endl;
}

const bool InterfaceTracking_Viewer::load_startframe() {
	frame = 0;
	return load(filename);
}

const bool InterfaceTracking_Viewer::load_endframe() {
	unsigned int endframe = 0;

	while(true) {
		string full_filename = get_concat_filename(filename, endframe);
		ifstream ifs(StringLibrary::string2wstring(full_filename.c_str()).c_str(), ios::binary);

		if(!ifs)
			break;

		endframe++;
	}

	frame = endframe - 1;
	return load(filename);
}

const bool InterfaceTracking_Viewer::load(const string& _filename) {
	filename = _filename;
	bool success = true;

	string full_filename = filename;

	ifstream ifs(StringLibrary::string2wstring(full_filename.c_str()).c_str(), ifstream::in | ios::binary);
	if(!ifs) {
		ifs.clear();

		full_filename = get_concat_filename(filename, frame);
		ifs.open(StringLibrary::string2wstring(full_filename.c_str()).c_str(), ifstream::in | ios::binary);

		if(!ifs) {
//			throwError("No file %s\n", full_filename.c_str());
			return false;
		}
	}

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(tracking == NULL) {
		if(InterfaceTracking_SRL2D::XmlElementName().compare(name) == 0)		tracking = new InterfaceTracking_SRL2D();
		else if(InterfaceTracking_SRL3D::XmlElementName().compare(name) == 0)	tracking = new InterfaceTracking_SRL3D();

		else if(InterfaceTracking_SLC3D::XmlElementName().compare(name) == 0)	tracking = new InterfaceTracking_SLC3D();

		else if(InterfaceTracking_Levelset2D::XmlElementName().compare(name) == 0) tracking = new InterfaceTracking_Levelset2D();
		else if(InterfaceTracking_Levelset3D::XmlElementName().compare(name) == 0) tracking = new InterfaceTracking_Levelset3D();

		else if(InterfaceTracking_PLS2D::XmlElementName().compare(name) == 0) tracking = new InterfaceTracking_PLS2D();
		else if(InterfaceTracking_PLS3D::XmlElementName().compare(name) == 0) tracking = new InterfaceTracking_PLS3D();
	}	

	if(tracking) {
		if(tracking->name_class().compare(name) != 0)
			return false;

		if(tracking->load(filename) == false) 
			success = tracking->load_compressed(get_concat_filename(filename, frame));
		
		if(success) {						
			log_info << "Liquid Volume : " << tracking->getTotalVolume() << endl;
			log_info << "InterfaceTracking_Viewer Load() Filename : " << filename.c_str() << " _ Frame : " << frame << endl;
		}		

		return success;
	}	
	return false;
}

const bool InterfaceTracking_Viewer::load_frame(int _frame) {
	const int pre_frame = frame;

	frame = _frame;

	if(load(filename) == false) {
		frame = pre_frame;
		return false;
	}

	return true;
}

const bool InterfaceTracking_Viewer::move_frame(int icr) {
	frame += icr;

	if(load(filename) == false) {
		cout << "load failed" << endl;

		frame -= icr;	
		return false;
	}

	return true;
}

void InterfaceTracking_Viewer::update(real dt) {
	isPause = !(move_frame(1));
	FluidSimulator::frame = frame;
}

void InterfaceTracking_Viewer::glBind() const {
	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	glDraw();	
}

void InterfaceTracking_Viewer::glDraw() const {

	if(tracking == NULL) return;

	Material mat;
	mat.setAmbient(RGBA::blue());
	mat.setDiffuse(RGBA::blue());
	mat.setSpecular(RGBA::blue() * 0.25f);
	mat.setShininess(8.5f);

	if(dLiquid) { 
		if(dMesh) {	tracking->glDrawInterface(mat, false);	}
		else	  {	tracking->glDrawInterface(mat, true);	}

		if(dInnerBdry) tracking->glDrawInnerBoundary();
	} 
	if(dPhi && tracking->getDimension() == 2 && tracking->isEulerian())	{ tracking->glDrawPhi();	}

	if(dAuxiliary) {
		InterfaceTracking_PLS2D* pls = dynamic_cast<InterfaceTracking_PLS2D*>(tracking);
		if(pls != NULL)	pls->glDrawParticles();		

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
const Vector3f InterfaceTracking_Viewer::viewEye() const { 
	if(tracking)	return tracking->viewEye();
	return Vector3f(0,0,1);
}

const Vector3f InterfaceTracking_Viewer::viewRef() const { 
	if(tracking)	return tracking->viewRef();
	return Vector3f(0,0,0);
}

const Vector3f InterfaceTracking_Viewer::viewUp() const  { 
	if(tracking)	return tracking->viewUp();
	return Vector3f(0,1,0);
}

const float InterfaceTracking_Viewer::viewTheta() const  { 
	if(tracking)	return tracking->viewTheta();
	return 45;
}

InterfaceTracking_Viewer_IOHandler::InterfaceTracking_Viewer_IOHandler(InterfaceTracking_Viewer *_s) {
	viewer = _s;
}

InterfaceTracking_Viewer_IOHandler::~InterfaceTracking_Viewer_IOHandler() {
}

void InterfaceTracking_Viewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'i' : viewer->dInnerBdry = !viewer->dInnerBdry; break;
		case 'm' : viewer->dMesh = !viewer->dMesh; break;
		case 'b' : viewer->dBoundingBox = !viewer->dBoundingBox; break;
		case 'g' : viewer->dGrid = !viewer->dGrid; break;
		case 'l' : viewer->dLiquid = !viewer->dLiquid; break;		
		case 'p' : viewer->dPhi = !viewer->dPhi; break;
		case 'a' : viewer->dAuxiliary = !viewer->dAuxiliary; break;	
		case 'x' : viewer->exportOBJ(get_concat_filename(output_folder + "interfaceTracking ", viewer->frame) + ".obj"); break;
		case 's' : viewer->load_startframe(); break;
		case 'e' : viewer->load_endframe(); break;
		case '*' : viewer->move_frame(+10); break;
		case '/' : viewer->move_frame(-10); break;
		case '+' : viewer->move_frame(+1); break;
		case '-' : viewer->move_frame(-1); break;		
		case ' ' : viewer->isPause = !viewer->isPause; break;

		case 'o' : scene->getCamera()->set(viewer->viewEye(), viewer->viewRef(), viewer->viewUp(), viewer->viewTheta()); break;
	}
}

void InterfaceTracking_Viewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void InterfaceTracking_Viewer_IOHandler::mouseMotionHandler(int x, int y) {
}
