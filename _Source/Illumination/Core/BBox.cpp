#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/BBox.h>

using namespace illuminationspace;

const bool BBox::rayIntersect(const Ray& r) const {
	real interval_min = r.tMin();
	real interval_max = r.tMax();

	int posneg = r.posneg[0];
	real t0 = (pp[posneg].x() - r.data[0].x()) * r.data[2].x();
	real t1 = (pp[1-posneg].x() - r.data[0].x()) * r.data[2].x();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;
	if(interval_min > interval_max) return false;

	posneg = r.posneg[1];
	t0 = (pp[posneg].y() - r.data[0].y()) * r.data[2].y();
	t1 = (pp[1-posneg].y() - r.data[0].y()) * r.data[2].y();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;
	if(interval_min > interval_max) return false;

	posneg = r.posneg[2];
	t0 = (pp[posneg].z() - r.data[0].z()) * r.data[2].z();
	t1 = (pp[1-posneg].z() - r.data[0].z()) * r.data[2].z();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;
	return (interval_min <= interval_max);
}

const bool BBox::intersectP(const Ray& r, real &_t0, real &_t1) const {
	real interval_min = r.tMin();
	real interval_max = r.tMax();

	int posneg = r.posneg[0];
	real t0 = (pp[posneg].x() - r.data[0].x()) * r.data[2].x();
	real t1 = (pp[1-posneg].x() - r.data[0].x()) * r.data[2].x();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;
	if(interval_min > interval_max) return false;

	posneg = r.posneg[1];
	t0 = (pp[posneg].y() - r.data[0].y()) * r.data[2].y();
	t1 = (pp[1-posneg].y() - r.data[0].y()) * r.data[2].y();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;
	if(interval_min > interval_max) return false;

	posneg = r.posneg[2];
	t0 = (pp[posneg].z() - r.data[0].z()) * r.data[2].z();
	t1 = (pp[1-posneg].z() - r.data[0].z()) * r.data[2].z();

	if(t0 > interval_min) interval_min = t0;
	if(t1 < interval_max) interval_max = t1;	

	_t0 = interval_min;
	_t1 = interval_max;
	return (interval_min <= interval_max);
}