//#ifndef _WIN64
#define _USE_XML_
//#endif

#ifdef _USE_XML_

#ifndef _BASIS_XML_LIBRARY_H_
#define _BASIS_XML_LIBRARY_H_

#include <iostream>
#include <string>
#include <libxml/parser.h>
#include "../Basis.h"
#include "../Basis/BasicLibrary.h"

namespace basis {

#define XML_GET_FUNCTION_HEADER(TYPE, FUNCNAME)											\
public :																				\
	static inline const TYPE FUNCNAME(xmlNodePtr node, const char* const keyWord) {		\
		return FUNCNAME(node, BAD_CAST keyWord);										\
	}																					\
private :																				\
	template<typename XMLPTR>															\
	static inline const TYPE FUNCNAME(XMLPTR node);										\
	static const TYPE FUNCNAME(xmlNodePtr node, const xmlChar* const keyWord);

#define XML_ADDCHILD_FUNCTION_HEADER(TYPE)																\
public :																								\
	static inline void xmlAddChild(xmlNodePtr parent, const char* const keyWord, const TYPE &data) {	\
		xmlNodePtr node = xmlNewNode(NULL, BAD_CAST keyWord);											\
		xmlPutProperty(node, "value", data);															\
		::xmlAddChild(parent, node);																	\
	}

#define XML_PUTPROP_FUNCTION_HEADER(TYPE, SPECIFIER)													\
public :																								\
	static inline void xmlPutProperty(xmlNodePtr node, const char* const keyWord, const TYPE &data) {	\
		static char _data[256];																			\
		sprintf_s(_data, SPECIFIER, data);																\
		xmlNewProp(node, BAD_CAST keyWord, BAD_CAST _data);												\
	}

#define XML_PUTTEXT_FUNCTION_HEADER(TYPE, SPECIFIER)												\
public :																							\
	static inline void xmlPutText(xmlNodePtr node, const TYPE &data) {								\
		static char _data[256];																		\
		sprintf_s(_data, SPECIFIER, data);															\
		xmlNodePtr textNode	= xmlNewText(BAD_CAST _data);											\
		::xmlAddChild(node, textNode);																\
	}																								\
																									\
	static inline void xmlPutText(xmlNodePtr node, const char* const keyWord, const TYPE &data) {	\
		xmlNodePtr newNode = xmlNewNode(NULL, BAD_CAST keyWord);									\
		xmlPutText(newNode, data);																	\
		::xmlAddChild(node, newNode);																\
	}
		
class BASISDLL XmlLibrary {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Import/Export XML File
	template<typename T> static void importXmlFile(T *desc, const char* const filename);
	template<typename T> static void exportXmlFile(const T *desc, const char* const filename);

	//////////////////////////////////////////////////////////////////////////
	// Find XML File/Attribute
	static inline const xmlNodePtr xmlFindNode(xmlDocPtr doc, const char* const keyWord);
	static inline const xmlNodePtr xmlFindNode(xmlNodePtr node, const char* const keyWord);
	static inline const xmlAttrPtr xmlFindAttr(xmlNodePtr node, const char* const keyWord);

public : 
	//////////////////////////////////////////////////////////////////////////
	// Get XML Value
	//static inline const TYPE xmlGet###(xmlNodePtr node, const char* const keyWord)
	XML_GET_FUNCTION_HEADER(bool, xmlGetBool);
	XML_GET_FUNCTION_HEADER(int, xmlGetInt);
	XML_GET_FUNCTION_HEADER(float, xmlGetFloat);
	XML_GET_FUNCTION_HEADER(double, xmlGetDouble);
	XML_GET_FUNCTION_HEADER(std::string, xmlGetString);

public : 
	//////////////////////////////////////////////////////////////////////////
	// Add XML Child
	template<typename T> static inline void xmlAddChild(xmlNodePtr parent, const char* const keyWord, const T& data);

	//static inline void xmlAddChild(xmlNodePtr parent, const char* const keyWord, const T &data);
	XML_ADDCHILD_FUNCTION_HEADER(bool);
	XML_ADDCHILD_FUNCTION_HEADER(int);
	XML_ADDCHILD_FUNCTION_HEADER(float);
	XML_ADDCHILD_FUNCTION_HEADER(double);
	XML_ADDCHILD_FUNCTION_HEADER(char*);	

public : 
	//////////////////////////////////////////////////////////////////////////
	// Put XML Properties
	static inline void xmlPutProperty(xmlNodePtr node, const char* const keyWord, const bool data);		
	XML_PUTPROP_FUNCTION_HEADER(int, "%d");
	XML_PUTPROP_FUNCTION_HEADER(float, "%f");
	XML_PUTPROP_FUNCTION_HEADER(double, "%g");
	XML_PUTPROP_FUNCTION_HEADER(char*, "%s");

public : 
	//////////////////////////////////////////////////////////////////////////
	// Put XML Text
	//static inline void xmlPutText(xmlNodePtr node, const bool data);		
	//static inline void xmlPutText(xmlNodePtr node, const char* const keyWord, const bool data);
	//XML_PUTTEXT_FUNCTION_HEADER(int, "%d");
	//XML_PUTTEXT_FUNCTION_HEADER(float, "%f");
	//XML_PUTTEXT_FUNCTION_HEADER(double, "%g");
	//XML_PUTTEXT_FUNCTION_HEADER(char*, "%s");

private :
	//////////////////////////////////////////////////////////////////////////
	// Find XML File/Attribute
	static const xmlNodePtr xmlFindNode(xmlDocPtr doc, const xmlChar* const keyWord);
	static const xmlNodePtr xmlFindNode(xmlNodePtr node, const xmlChar* const keyWord);
	static const xmlNodePtr xmlFindNodeRecursion(xmlNodePtr node, const xmlChar* const keyWord);	
	static const xmlAttrPtr xmlFindAttr(xmlNodePtr node, const xmlChar* const keyWord);		

	//////////////////////////////////////////////////////////////////////////
	// Add XML Child
	template<typename T> static inline void xmlAddChild(xmlNodePtr parent, const xmlChar* const keyWord, const T& data);
};

//////////////////////////////////////////////////////////////////////////
// Import/Export XML File
template<typename T> 
void XmlLibrary::importXmlFile(T *desc, const char* const filename) {

	xmlDocPtr	doc		= NULL;
	xmlNodePtr	node	= NULL;

	// Read file
	doc = xmlReadFile(filename, NULL, 0);
	if(doc == NULL) {
		throwError("XML file open error : %s\n", filename);
		exit(1);
	}
	if(desc == NULL)	{
		throwError("Class is not constructed yet : %s\n", T::XmlElementName().c_str());
		exit(1);
	}
	if((node = xmlFindNode(doc, T::XmlElementName().c_str())) == NULL) {
		throwError("Cannot Find XmlNode %s in File %s\n", T::XmlElementName().c_str(), filename);
		exit(1);
	}	

	desc->importXmlNode(node);

	// Memory reaping
	xmlFreeDoc(doc);
	xmlCleanupParser();
}

template<typename T> 
void XmlLibrary::exportXmlFile(const T *desc, const char* const filename) {

	xmlDocPtr	doc		= NULL;
	xmlNodePtr	node	= NULL;

	// Root node
	doc = xmlNewDoc(BAD_CAST "1.0");
	if((node = desc->exportXmlNode()) == NULL){
		throwError("Cannot Export XmlNode %s\n", T::XmlElementName().c_str());
		exit(1);
	}

	xmlDocSetRootElement(doc, node);	

	// Save file
	xmlSaveFormatFileEnc(filename, doc, "UTF-8", 1);

	// Memory reaping
	xmlFreeDoc(doc);
	xmlCleanupParser();
	xmlMemoryDump();	
}

//////////////////////////////////////////////////////////////////////////
// Find XML File/Attribute
inline const xmlNodePtr XmlLibrary::xmlFindNode(xmlDocPtr doc, const char* const keyWord) {
	return xmlFindNode(doc, BAD_CAST keyWord);
}

inline const xmlNodePtr XmlLibrary::xmlFindNode(xmlNodePtr node, const char* const keyWord) {
	return xmlFindNode(node, BAD_CAST keyWord);
}

inline const xmlAttrPtr XmlLibrary::xmlFindAttr(xmlNodePtr node, const char* const keyWord) {
	return xmlFindAttr(node, BAD_CAST keyWord);
}

//////////////////////////////////////////////////////////////////////////
// Add XML Child
template<typename T> 
inline void XmlLibrary::xmlAddChild(xmlNodePtr parent, const char* const keyWord, const T& data) {
	xmlAddChild(parent, BAD_CAST keyWord, data);
}

template<typename T> 
inline void XmlLibrary::xmlAddChild(xmlNodePtr parent, const xmlChar* const keyWord, const T& data) {

	xmlNodePtr node = xmlNewNode(NULL, keyWord);
	::xmlAddChild(node, data.exportXmlNode());
	::xmlAddChild(parent, node);
}

//////////////////////////////////////////////////////////////////////////
// Get XML Value
template<typename XMLPTR> 
inline const bool XmlLibrary::xmlGetBool(XMLPTR node) { 
	if(xmlStrcmp(node->children->content, BAD_CAST "true")==0) 
		return true;
	else
		return false;	
}

template<typename XMLPTR> 
inline const int XmlLibrary::xmlGetInt(XMLPTR node) { 
	return atoi((char*)node->children->content);		
}

template<typename XMLPTR> 
inline const float XmlLibrary::xmlGetFloat(XMLPTR node)	{ 
	return (float)atof((char*)node->children->content);	
}

template<typename XMLPTR> 
inline const double XmlLibrary::xmlGetDouble(XMLPTR node)	{ 
	return atof((char*)node->children->content);	
}

template<typename XMLPTR> 
inline const std::string XmlLibrary::xmlGetString(XMLPTR node) { 
	return std::string((char*)node->children->content);		
}

//////////////////////////////////////////////////////////////////////////
// Put XML Properties
inline void XmlLibrary::xmlPutProperty(xmlNodePtr node, const char* const keyWord, const bool data) {	
	if(data == true)	xmlNewProp(node, BAD_CAST keyWord, BAD_CAST "true");
	else				xmlNewProp(node, BAD_CAST keyWord, BAD_CAST "false");
}

//////////////////////////////////////////////////////////////////////////
// Put XML Text
//inline void XmlLibrary::xmlPutText(xmlNodePtr node, const bool data) {
//	if(data == true)	::xmlAddChild(node, xmlNewText(BAD_CAST "true"));
//	else				::xmlAddChild(node, xmlNewText(BAD_CAST "false"));			
//}	
//inline void XmlLibrary::xmlPutText(xmlNodePtr node, const char* const keyWord, const bool data) {	
//	xmlNodePtr newNode = xmlNewNode(NULL, BAD_CAST keyWord);								
//	xmlPutText(newNode, data);																
//	::xmlAddChild(node, newNode);																
//}

#undef XML_GET_FUNCTION_HEADER
#undef XML_ADDCHILD_FUNCTION_HEADER
#undef XML_PUTPROP_FUNCTION_HEADER
#undef XML_PUTTEXT_FUNCTION_HEADER

};		// namespace basis

#endif

#endif
