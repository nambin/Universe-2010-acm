#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/SimpleAggregate.h>

using namespace illuminationspace;

const BBox SimpleAggregate::SimpleAggregateValue::worldBound() const {

	BBox bbox;
	for(unsigned int i=0;i<prims.size();++i) 
		bbox.surround(prims[i].worldBound());

	return bbox;
}

const bool SimpleAggregate::SimpleAggregateValue::hit(const Ray& r, Intersection &intersection) const {

	bool is_a_hit = false;

	for(unsigned int i=0;i<prims.size();++i) {
		if(prims[i].hit(r, intersection)) {
			r.tmax = intersection.tHit;
			is_a_hit = true;
		}
	}

	return is_a_hit;
}

const bool SimpleAggregate::SimpleAggregateValue::shadowHit(const Ray& r) const {

	bool is_a_hit = false;

	for(unsigned int i=0;i<prims.size();++i) {
		if(prims[i].shadowHit(r))
			return true;
	}

	return is_a_hit;
}

void SimpleAggregate::SimpleAggregateValue::addPrimitive(Primitive prim) {
	//prims.push_back(prim);

	std::vector<Primitive> _prims;	
	prim.fullyRefine(_prims);

	for(unsigned int i=0;i<_prims.size();++i) 
		prims.push_back(_prims[i]);
}

void SimpleAggregate::SimpleAggregateValue::addPrimitive(const std::vector<Primitive> &_prims) {
	for(unsigned int i=0;i<_prims.size();++i)
		addPrimitive(_prims[i]);
}