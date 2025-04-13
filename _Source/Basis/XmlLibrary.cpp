#include <stdafx.h>
#include <Basis/XmlLibrary/XmlLibrary.h>

using namespace std;
using namespace basis;

#ifdef _USE_XML_

//////////////////////////////////////////////////////////////////////////
// Find XML File/Attribute
const xmlNodePtr XmlLibrary::xmlFindNode(xmlDocPtr doc, const xmlChar* const keyWord){
	xmlNodePtr rootNode = xmlDocGetRootElement(doc);
	return xmlFindNodeRecursion(rootNode, keyWord);
}

const xmlNodePtr XmlLibrary::xmlFindNode(xmlNodePtr node, const xmlChar* const keyWord){
	if(node == NULL)
		return NULL;

	xmlNodePtr childNode = NULL;

	if(node->type == XML_ELEMENT_NODE)
		if(xmlStrcmp(node->name, keyWord)==0)
			return node;

	childNode = xmlFindNodeRecursion(node->children, keyWord);
	if(childNode)
		return childNode;

	return NULL;
}

const xmlNodePtr XmlLibrary::xmlFindNodeRecursion(xmlNodePtr node, const xmlChar* const keyWord) {
	xmlNodePtr curNode		= NULL;
	xmlNodePtr childNode	= NULL;

	for(curNode = node; curNode; curNode = curNode->next) 
		if(curNode->type == XML_ELEMENT_NODE) 
			if(xmlStrcmp(curNode->name, keyWord)==0) 
				return curNode;

	for(curNode = node; curNode; curNode = curNode->next) {
		childNode = xmlFindNodeRecursion(curNode->children, keyWord);
		if(childNode) 
			return childNode;
	}

	return NULL;
}

const xmlAttrPtr XmlLibrary::xmlFindAttr(xmlNodePtr node, const xmlChar* const keyWord) {
	if(node == NULL)
		return NULL;

	xmlAttrPtr rtnAttrNode;	
	for(rtnAttrNode = node->properties; rtnAttrNode; rtnAttrNode = rtnAttrNode->next) 
		if(xmlStrcmp(rtnAttrNode->name, keyWord)==0) 
			return rtnAttrNode;

	return NULL;
}

//////////////////////////////////////////////////////////////////////////
// Get XML Value
#define XML_GET_FUNCTION_SOURCE(TYPE, FUNCNAME, ERRORRET)							\
const TYPE XmlLibrary::FUNCNAME(xmlNodePtr node, const xmlChar * const keyWord) {	\
																					\
	xmlAttrPtr attr = xmlFindAttr(node, keyWord);									\
	if(attr != NULL) return FUNCNAME(attr);											\
																					\
	xmlNodePtr curNode = xmlFindNode(node, keyWord);								\
	if(curNode != NULL) {															\
		attr = xmlFindAttr(curNode, "value");										\
		if(attr != NULL)															\
			return FUNCNAME(attr);													\
																					\
		return FUNCNAME(curNode);													\
	}																				\
	else {																			\
		throwError("XML Get Property error : %s", (char*)keyWord);					\
		return ERRORRET;															\
	}																				\
}

XML_GET_FUNCTION_SOURCE(bool, xmlGetBool, false)
XML_GET_FUNCTION_SOURCE(int, xmlGetInt, 0)
XML_GET_FUNCTION_SOURCE(float, xmlGetFloat, 0.0f)
XML_GET_FUNCTION_SOURCE(double, xmlGetDouble, 0.0)
XML_GET_FUNCTION_SOURCE(string, xmlGetString, string(""));

#undef XML_GET_FUNCTION_SOURCE

#endif