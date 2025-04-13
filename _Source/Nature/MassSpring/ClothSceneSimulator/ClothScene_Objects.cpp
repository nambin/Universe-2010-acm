#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothScene_Objects.h>

using namespace springspace;

void ClothScene_Objects::push_DynamicObject(DynamicObject<Vector3<real>> obj) {
	objects.push_back(obj);
}

void ClothScene_Objects::query(QueryInfo& info, const Vector3<real>& pos) const {
	size_t idx = -1;
	info.phi = 1.0e+6;
	for(size_t i=0;i<objects.size();++i) {
		real phi = objects[i].getPhi(pos);
		if(phi < info.phi) {
			idx = i;
			info.phi = phi;
		}
	}

	assert(idx != -1);
	info.normal = objects[idx].getNormal(pos);
	info.vel = objects[idx].getVel(pos);
	info.mu = objects[idx].getFrictionCoefficient(pos);
}

void ClothScene_Objects::glDraw() const {
	for(size_t i=0;i<objects.size();++i)
		objects[i].glDraw();
}