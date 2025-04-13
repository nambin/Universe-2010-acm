#ifndef _ILLUMINATION_AGGREGATE_H_
#define _ILLUMINATION_AGGREGATE_H_

#include "Primitive.h"

namespace illuminationspace {

class ILLUMINATIONDLL Aggregate : public Primitive {
protected :
	class AggregateValue;

public :
	Aggregate() {}
	Aggregate(AggregateValue *ptr) : Primitive(ptr) {}

protected :
	class ILLUMINATIONDLL AggregateValue : public Primitive::PrimitiveValue {
	public : 
		AggregateValue() {}		
	};
};

struct MailboxPrim {
public :
	MailboxPrim(const Primitive &p)
	: primitive(p), lastMailboxId(-1) {}

public :
	Primitive primitive;
	int lastMailboxId;
};

};		// illuminationspace

#endif
