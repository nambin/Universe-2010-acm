#ifndef _ILLUMINATION_SIMPLE_AGGREGATE_H_
#define _ILLUMINATION_SIMPLE_AGGREGATE_H_

#include "Aggregate.h"

namespace illuminationspace {

class ILLUMINATIONDLL SimpleAggregate : public Aggregate {
protected : 
	class SimpleAggregateValue;

public : 
	SimpleAggregate() : Aggregate(new SimpleAggregate::SimpleAggregateValue()) {}
	SimpleAggregate(const std::vector<Primitive> &p) : Aggregate(new SimpleAggregate::SimpleAggregateValue(p)) {}

	void addPrimitive(Primitive p) { 
		dynamic_cast<SimpleAggregateValue*>(const_cast<PrimitiveValue*>(value.getptr()))->addPrimitive(p); 
	}
	void addPrimitive(const std::vector<Primitive> &prims)	{ 
		dynamic_cast<SimpleAggregateValue*>(const_cast<PrimitiveValue*>(value.getptr()))->addPrimitive(prims); 
	}

protected :
	class ILLUMINATIONDLL SimpleAggregateValue : public Aggregate::AggregateValue {
	public : 
		SimpleAggregateValue() {}
		SimpleAggregateValue(const std::vector<Primitive> &p) { addPrimitive(p); }		

		SimpleAggregateValue * clone() const { return new SimpleAggregateValue(*this); }

		const BBox worldBound() const;
		const bool canIntersect() const { return false; }

		const bool hit(const Ray& r, Intersection &intersection) const;
		const bool shadowHit(const Ray& r) const;

		void addPrimitive(Primitive);
		void addPrimitive(const std::vector<Primitive> &prims);

	public : 
		std::vector<Primitive> prims;
	};	
};

};		// illuminationspace

#endif