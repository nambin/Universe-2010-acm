#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/AreaLight.h>
#include <Illumination/IlluminationHeader/Illumination/Random.h>
#include <Illumination/IlluminationHeader/Core/IlluminationLibrary.h>

using namespace illuminationspace;

AreaLight::AreaLightValue::AreaLightValue(const Transform &l2w, const Spectrum &le, const Shape s, int ns, int ns_stratified, Sampler::SAMPLING_TYPE type)	
: Light::LightValue(l2w, ns, ns_stratified, type) {
	_L = le;
	if(s.canIntersect())	
		shape = s;
	else {
		std::vector<Shape> todo, done;
		todo.push_back(s);
		while (todo.size()) {
			Shape sh = todo.back();
			todo.pop_back();
			if (sh.canIntersect())	done.push_back(sh);
			else					sh.refine(todo);
		}
		if (done.size() == 1) shape = done[0];
		else {
			if (done.size() > 16)
				std::cout << "Area light geometry turned into %d shapes; " << done.size() << " It may be very inefficient." << std::endl;
			shape = ShapeSet(done, s.getShapeValue()->ObjectToWorld);
		}
	}

	_area = shape.area();
}

const real AreaLight::AreaLightValue::pdf(const Vector3<real> &pt, const Vector3<real> &wi) const	{
	return shape.pdf();
}

const Spectrum AreaLight::AreaLightValue::sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const	{
	Vector3<real> ns;
	Vector3<real> ps = shape.sample(0.5f, 0.5f, &ns);

	*wi = unitVector(ps - pt);
	vis->setSegment(pt, ps);
	return L(pt, ns, -(*wi)) / (shape.pdf() * sqrlength(ps, pt)) * abs(dot(ns, *wi));
}

const Spectrum AreaLight::AreaLightValue::sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const {
	Vector3<real> ns;
	Vector3<real> ps = shape.sample(u1, u2, &ns);

	*wi = unitVector(ps - pt);
	*pdf = shape.pdf() * sqrlength(ps, pt) / abs(dot(ns, *wi));
	vis->setSegment(pt, ps);
	return L(pt, ns, -(*wi));
}

void AreaLight::AreaLightValue::generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const {
	const ONB &uvw = shape.sample(random(), random(), pos);
	Vector3<real> wi = IlluminationLibrary::CosineSampleHemisphere(random(), random());	

	dir = wi[0] * uvw.u() + wi[1] * uvw.v() + wi[2] * uvw.w();	
}
