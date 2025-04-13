#include <stdafx.h>
#include <vector>
#include <Illumination/IlluminationHeader/PhotonMap/PhotonEmitter.h>

using namespace std;
using namespace illuminationspace;

PhotonEmitter_UnbiasedEmitter::PhotonEmitter_UnbiasedEmitter(size_t num, PhotonMap* global, PhotonMap* caustic) 
: numPhotons(num), globalPhotonMap(global), causticPhotonMap(caustic)
{
}

void PhotonEmitter_UnbiasedEmitter::emitPhotons(const Scene* scene) const {
	real totalPower = 0;
	for(size_t i=0;i<scene->lights.size();++i)	
		totalPower += scene->lights[i].power(scene).sum();

	vector<size_t> nPhotonsPerLight;
	for(size_t i=0;i<scene->lights.size();++i)
		nPhotonsPerLight.push_back((scene->lights[i].power(scene).sum() / totalPower * numPhotons));

	Vector3<real> pos, dir;
	for(size_t i=0;i<scene->lights.size();++i) {

		size_t nTotalPhotons = 0, nValidPhotons = 0;
		while(nValidPhotons < nPhotonsPerLight[i])	{

			scene->lights[i].generatePhotonPath(pos, dir);

			Ray ray(pos, dir, Spectrum(0.0f), 0);
			Spectrum power = scene->lights[i].power(scene);
			if(tracePhoton(scene, ray, power, globalPhotonMap, causticPhotonMap))
				nValidPhotons++;
			nTotalPhotons++;			
		}

		globalPhotonMap->scale_photon_power(1.0f/nTotalPhotons);
		causticPhotonMap->scale_photon_power(1.0f/nTotalPhotons);
	}

	globalPhotonMap->balance();
	causticPhotonMap->balance();

	cout << "Total Global  Photons : " << globalPhotonMap->nStoredPhotons() << endl;
	cout << "Total Caustic Photons : " << causticPhotonMap->nStoredPhotons() << endl;
}

bool PhotonEmitter_UnbiasedEmitter::tracePhoton(const Scene* scene, Ray& r, Spectrum& power, PhotonMap *globalPhotonMap, PhotonMap *causticPhotonMap) {
	Intersection intersection;	
	Spectrum retPower;
	BxDFType reflectType;
	Vector3<real> wi;
	real _power[3], initialPower3 = power.sum();
	int count = 0;

	bool isHitted = false;
	bool isHitDiffuse = false, isHitSpecular = false;
	bool isHitSpecular_JustBefore = false;

	while(true) {

		if(!scene->hit(r, intersection))	
			return (isHitted || false);
		isHitted = true;

		BSDF *bsdf = intersection.getBSDF(r);

		if(bsdf == NULL)
			return true;

		if(bsdf->numComponents(BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_TRANSMISSION)) != 0) {
			_power[0] = power[0];
			_power[1] = power[1];
			_power[2] = power[2];

			if(!isHitDiffuse && isHitSpecular)
				//if(isHitSpecular_JustBefore)
				causticPhotonMap->store(_power, intersection.dg.hit_p.data, r.direction().data);

			globalPhotonMap->store(_power, intersection.dg.hit_p.data, r.direction().data);
		}

		if(!bsdf->reflectPhoton(r.direction(), &wi, random(), random(), random(), power, &retPower, &reflectType)) {
			BSDF::FreeAll();
			return true;
		}

		if(count > 20 || retPower.sum() < 1.0e-2 * initialPower3)
			return true;

		BSDF::FreeAll();

		isHitDiffuse = (isHitDiffuse || (reflectType & BSDF_DIFFUSE));
		isHitSpecular = (isHitSpecular || (reflectType & BSDF_SPECULAR));		
		isHitSpecular_JustBefore = (true && (reflectType & BSDF_SPECULAR));

		r.setRay(intersection.dg.hit_p, wi);
		power = retPower;		

		count++;
	}		
}