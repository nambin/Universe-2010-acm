#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothScene.h>

using namespace std;
using namespace springspace;

void ClothScene::glBind() const {
	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDraw();	
}

void ClothScene::glDraw() const {	
	if(dSurface)		mesh->glDraw_SmoothSurface(ClothLibrary::Material_Front_Basic(), ClothLibrary::Material_Back_Basic());
	if(dParticle)		mesh->glDraw_Particle(RGBA::red(), RGBA::yellow(), mesh->width() * 0.01f);
	if(dStretchSpring)	mesh->glDraw_StretchSpring(RGBA::white(), 2.0f);
	if(dShearSpring)	mesh->glDraw_ShearSpring(RGBA::white(), 2.0f);
	if(dBendingSpring)	mesh->glDraw_BendingSpring(RGBA::white(), 2.0f);

	if(dObjects && objects)	objects->glDraw();		
}

//////////////////////////////////////////////////////////////////////////
// Initial View Setting
const Vector3f ClothScene::viewEye() const { return Vector3f(0, mesh->width() * 0.3f, MAX(mesh->width(), mesh->length()) * 2.5f);	}
const Vector3f ClothScene::viewRef() const { return Vector3f(0, 0, 0);												}
const Vector3f ClothScene::viewUp() const  { return Vector3f(0, 1, 0);												}
const float ClothScene::viewTheta() const  { return 45.0f;															}
