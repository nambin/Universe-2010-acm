#include <stdafx.h>
#include <Nature/Fluid/FluidScene/FluidScene_Objects.h>

using namespace nature;
using namespace fluidspace;

void FluidScene_Objects2D::push_DynamicObject(DynamicObject<Vector2<real>> obj) {
	objects.push_back(obj);
}

const bool FluidScene_Objects2D::isInside(const Vector2<real>& p) const {
	for(size_t i=0;i<objects.size();++i) {
		if(objects[i].isInside(p))
			return true;
	}
	return false;
}

const Vector2<real> FluidScene_Objects2D::getVel(const Vector2<real>& p) const {
	for(size_t i=0;i<objects.size();++i) {
		if(objects[i].isInside(p))
			return objects[i].getVel(p);
	}

	real min_phi = 1.0e+8; size_t idx;

	for(size_t i=0;i<objects.size();++i) {
		real phi = objects[i].getPhi(p);
		if(phi < min_phi) {
			idx = i; 
			min_phi = phi;
		}			
	}

	return objects[idx].getVel(p);
}

void FluidScene_Objects2D::glDraw() const {
	if(dObject == false)
		return;

	for(size_t i=0;i<objects.size();++i)
		objects[i].glDraw();
}