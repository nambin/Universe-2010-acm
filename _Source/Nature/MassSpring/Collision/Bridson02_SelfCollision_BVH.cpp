#include <stdafx.h>
#include <Nature/MassSpring/Collision/Bridson02_SelfCollision.h>

using namespace std;
using namespace geometryspace;
using namespace springspace;

BVH<unsigned int, real>* Bridson02_SelfCollision::makeStaticBVH(const Vec3Array<real>& pos, const Array<pair<unsigned int, unsigned int>>& order, const vector<unsigned int>& v2f, real eps) {
	unsigned int numFaces = static_cast<unsigned int>(v2f.size() / 3);

	Array<BVH<unsigned int, real>*> bvhPtr(numFaces);
	for(unsigned int i=0;i<numFaces;++i) {
		AABB<real> bbox;
		bbox.surround(pos[v2f[i*3+0]]);
		bbox.surround(pos[v2f[i*3+1]]);
		bbox.surround(pos[v2f[i*3+2]]);
		bbox.expand(eps);

		bvhPtr[i] = new BVH<unsigned int, real>(bbox);
		bvhPtr[i]->pushData(i);
	}

	return mergeBVH(order, bvhPtr);
}

BVH<unsigned int, real>* Bridson02_SelfCollision::makeDynamicBVH(const Vec3Array<real>& pos, const Vec3Array<real>& vel, const Array<pair<unsigned int, unsigned int>>& order, const vector<unsigned int>& v2f, real eps) {
	unsigned int numFaces = static_cast<unsigned int>(v2f.size() / 3);

	Array<BVH<unsigned int, real>*> bvhPtr(numFaces);
	for(unsigned int i=0;i<numFaces;++i) {
		AABB<real> bbox;
		bbox.surround(pos[v2f[i*3+0]]);
		bbox.surround(pos[v2f[i*3+1]]);
		bbox.surround(pos[v2f[i*3+2]]);
		bbox.surround(pos[v2f[i*3+0]] + vel[v2f[i*3+0]]);
		bbox.surround(pos[v2f[i*3+1]] + vel[v2f[i*3+1]]);
		bbox.surround(pos[v2f[i*3+2]] + vel[v2f[i*3+2]]);
		bbox.expand(eps);

		bvhPtr[i] = new BVH<unsigned int, real>(bbox);
		bvhPtr[i]->pushData(i);
	}

	return mergeBVH(order, bvhPtr);
}

BVH<unsigned int, real>* Bridson02_SelfCollision::mergeBVH(const Array<pair<unsigned int, unsigned int>>& order, Array<BVH<unsigned int, real>*>& bvhPtr) {
	size_t numBVH = bvhPtr.size();
	bool rightOrder = true;

	while(numBVH != 1) {
		for(size_t i=0;i<order.size();++i) {
			size_t idx = i;
			if(!rightOrder)
				idx = order.size() - 1 - i;

			if(bvhPtr[order[idx].first] == bvhPtr[order[idx].second])
				continue;

			BVH<unsigned int, real>* bvh = new BVH<unsigned int, real>(bvhPtr[order[idx].first], bvhPtr[order[idx].second], true);
			vector<unsigned int>& data = bvh->getData();

			for(size_t k=0;k<data.size();++k) 
				bvhPtr[data[k]] = bvh;

			--numBVH;
		}
		rightOrder = !rightOrder;
	}

	bvhPtr[0]->removeDuplicateData();
	return bvhPtr[0];
}
