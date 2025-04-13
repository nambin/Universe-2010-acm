#include <stdafx.h>
#include <Basis/XmlLibrary/XmlLibrary.h>
#include <Illumination/IlluminationHeader/Core/Transform.h>

using namespace std;
using namespace illuminationspace;

#ifdef _USE_XML_

void Transform::importXmlNode(xmlNodePtr node) {
	identity();
	Transform tmp;

	for(xmlNodePtr curNode = node->children; curNode; curNode = curNode->next){
		if(curNode->type == XML_ELEMENT_NODE) {
			tmp.getTransformFromXmlNode(curNode);
			(*this) = tmp * (*this);
		}			
	}	
}

void Transform::getTransformFromXmlNode(xmlNodePtr node) {
	xmlNodePtr curNode = NULL;
	if(xmlStrcmp(node->name, BAD_CAST "Translate")==0) {
		Vector3<real> v;
		v[0] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "x"));	
		v[1] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "y"));	
		v[2] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "z"));	
		translate(v);
		return;
	}
	if(xmlStrcmp(node->name, BAD_CAST "Rotate")==0) {
		Vector3<real> v; real angle;
		v[0] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "x"));	
		v[1] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "y"));	
		v[2] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "z"));	
		angle = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "angle"));	
		rotate(angle, v);
		return;
	}
	if(xmlStrcmp(node->name, BAD_CAST "Scale")==0) {
		Vector3<real> v;
		v[0] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "x"));	
		v[1] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "y"));	
		v[2] = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "z"));	
		scale(v[0], v[1], v[2]);
		return;
	}
}

#endif

const Ray Transform::transformRay(const Ray& r) const	{
	Ray rRet(r);

	rRet.setOrigin(transformPoint(r.data[0]));
	rRet.setDirection(transformVector(r.data[1]));
	return rRet;
}

void Transform::transformRay(const Ray& r, Ray* rTrans) const	{

	if(rTrans == NULL) {
		throwError("Error in Transform::tranformRay\n");
		exit(0);
	}

	rTrans->setRay(transformPoint(r.data[0]), transformVector(r.data[1]), r.extinction(), r.depth);
	rTrans->setTmin(r.tMin());
	rTrans->setTmax(r.tMax());
}

const BBox Transform::transformBBox(const BBox &bbox) const	{
	BBox ret    ( transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[0].y(), bbox.pp[0].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[0].y(), bbox.pp[1].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[1].y(), bbox.pp[0].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[0].y(), bbox.pp[0].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[1].y(), bbox.pp[0].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[0].y(), bbox.pp[1].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[1].y(), bbox.pp[1].z())) );
	ret.surround( transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[1].y(), bbox.pp[1].z())) );
	return ret;
}

const ONB Transform::transformONB(const ONB &uvw) const {
	ONB _uvw;
	_uvw.U = transformVector(uvw.U);
	_uvw.V = transformVector(uvw.V);
	_uvw.W = transformNormal(uvw.W);
	return _uvw;
}

//inline void Transform::transformBBox(const BBox &bbox, BBox &bboxTrans) const	{
//	bboxTrans.set( transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[0].y(), bbox.pp[0].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[0].y(), bbox.pp[1].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[1].y(), bbox.pp[0].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[0].y(), bbox.pp[0].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[1].y(), bbox.pp[0].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[0].y(), bbox.pp[1].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[0].x(), bbox.pp[1].y(), bbox.pp[1].z())) );
//	bboxTrans.surround(transformPoint(Vector3<real>(bbox.pp[1].x(), bbox.pp[1].y(), bbox.pp[1].z())) );
//}
