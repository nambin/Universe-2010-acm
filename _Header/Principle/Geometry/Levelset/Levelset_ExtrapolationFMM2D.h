#pragma once

namespace geometryspace {

template<typename T>
void Levelset_ExtrapolationFMM<T>::extrapolate(Grid2<T,T>& val, typename MacGrid2D<T>::FIELD_TYPE _type, const GridBoundaryCondition<T,T>& bdry, const Grid2<T,T>& _phi) const {
	Grid2<T,T> phi; phi.resize(val);
	Levelset2D<T>::getPhi(phi, _phi, _type);
	extrapolate(val, bdry, phi);
}

template<typename T>
void Levelset_ExtrapolationFMM<T>::extrapolate(Grid2<T,T>& val, const GridBoundaryCondition<T,T>& bdry, const Grid2<T,T>& phi) const {

	T h = MAX(phi.dx, phi.dy);
	grid_index_type _narrowBandWidth = MIN(static_cast<real>(MAX(phi.nX, phi.nY)), narrowBandWidth);
	minHeapSize = static_cast<grid_index_type>(3.0f * _narrowBandWidth * MAX(phi.nX,phi.nY));

//	if(heap == NULL)	heap = new FMMHeap<T>(minHeapSize, val.size());
//	else				heap->init(minHeapSize);

	if(heap) 
		delete heap;

	heap = new FMMHeap<T>(minHeapSize, val.size());

	Grid2<char,T> valTag(val.nX, val.nY, val.dx, val.dy, val.offset);

	for(grid_index_type i=0;i<val.nX;++i)	{
		for(grid_index_type j=0;j<val.nY;++j)	{
			if(phi(i,j) < 0)	valTag(i,j) = FMMHeap<T>::FMM_ACCEPTED;
			else				valTag(i,j) = FMMHeap<T>::FMM_FAR;
		}
	}

	for(grid_index_type i=0;i<val.nX;++i)	{
		for(grid_index_type j=0;j<val.nY;++j)	{
			if(valTag(i,j)==FMMHeap<T>::FMM_ACCEPTED)	{
				if(i != 0 && valTag(i-1,j) == FMMHeap<T>::FMM_FAR)	{
					valTag(i-1,j) = FMMHeap<T>::FMM_TRIAL;
					heap->push(val.idx(i-1,j), phi(i-1,j));
				}
				if(i != val.nX-1 && valTag(i+1,j) == FMMHeap<T>::FMM_FAR)	{
					valTag(i+1,j) = FMMHeap<T>::FMM_TRIAL;
					heap->push(val.idx(i+1,j), phi(i+1,j));
				}
				if(j != 0 && valTag(i,j-1) == FMMHeap<T>::FMM_FAR)	{
					valTag(i,j-1) = FMMHeap<T>::FMM_TRIAL;
					heap->push(val.idx(i,j-1), phi(i,j-1));
				}
				if(j != val.nY-1 && valTag(i,j+1) == FMMHeap<T>::FMM_FAR)	{
					valTag(i,j+1) = FMMHeap<T>::FMM_TRIAL;
					heap->push(val.idx(i,j+1), phi(i,j+1));
				}
			}
		}
	}

	while(!heap->empty())	{
		NodeHeap<T> nodeHeap = heap->pop();

		grid_index_type i = nodeHeap.idx / valTag.nY;
		grid_index_type j = nodeHeap.idx % valTag.nY;

		if(phi(i,j) > narrowBandWidth * h)
			break;

		valTag(nodeHeap.idx) = FMMHeap<T>::FMM_ACCEPTED;		

		real s,t, val_s, val_t;
		s = t = FMMHeap<T>::MAX_DISTANCE;

		if(i != 0 && valTag(i-1,j) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(phi(i-1,j) < s)	{
				val_s = val(i-1,j);
				s = phi(i-1,j);
			}
		}
		if(i != val.nX-1 && valTag(i+1,j) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(phi(i+1,j) < s)	{
				val_s = val(i+1,j);
				s = phi(i+1,j);
			}
		}
		if(j != 0 && valTag(i,j-1) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(phi(i,j-1) < t)	{
				val_t = val(i,j-1);
				t = phi(i,j-1);
			}
		}
		if(j != val.nY-1 && valTag(i,j+1) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(phi(i,j+1) < t)	{
				val_t = val(i,j+1);
				t = phi(i,j+1);
			}
		}

		if(s==FMMHeap<T>::MAX_DISTANCE)			val(nodeHeap.idx) = val_t;
		else if(t==FMMHeap<T>::MAX_DISTANCE)	val(nodeHeap.idx) = val_s;
		else {
			real ws = abs(phi(i,j)-s);
			real wt = abs(phi(i,j)-t);

			if(ws+wt < 1.e-5)	val(nodeHeap.idx) = 0.5f * (val_s+val_t);
			else				val(nodeHeap.idx) = (ws*val_s + wt*val_t) / (ws+wt);
		}

		if(i != 0 && valTag(i-1,j) == FMMHeap<T>::FMM_FAR)	{
			valTag(i-1,j) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i-1,j), phi(i-1,j));
		}
		if(i != val.nX-1 && valTag(i+1,j) == FMMHeap<T>::FMM_FAR)	{
			valTag(i+1,j) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i+1,j), phi(i+1,j));
		}
		if(j != 0 && valTag(i,j-1) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j-1) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j-1), phi(i,j-1));
		}
		if(j != val.nY-1 && valTag(i,j+1) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j+1) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j+1), phi(i,j+1));
		}
	}

	for(grid_index_type i=0;i<val.size();i++) {
		if(valTag(i) != FMMHeap<T>::FMM_ACCEPTED)	val(i) = 0.0f;
	}

	bdry.enforce(val);
}

};			// geometryspace