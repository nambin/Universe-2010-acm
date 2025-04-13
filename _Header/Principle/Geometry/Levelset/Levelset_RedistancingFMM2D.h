#pragma once

namespace geometryspace {

template<typename T>
void Levelset_RedistancingFMM<T>::redistance(Levelset2D<T>& lv) const {
	// Fast Marching in New NarrowBand
//	lv.phi.setNeumannBoundary();
	lv.phi.scale(lv.hInv);

	// Data Structure & 
	// Distance, sTag2D Initialization near interface
	redistance_initialize(lv);
//	return;

	// FMM initialization for positive, negative Region
	for(int k=0;k<2;k++)	{
		if(k==1)	{
			heap->init(minHeapSize);
		}

		for(size_t l=0;l<interfaceCellList.size();l+=2)	{
			grid_index_type i = interfaceCellList[l];
			grid_index_type j = interfaceCellList[l+1];

			if( lv.phi(i,j) < 0 || (lv.phi(i,j) == 0 && k==0))
				continue;		

			grid_index_type nearS;
			T nearPhi = 0;			

			if( i != 0 && sTag2D[(nearS = lv.iS(i-1,j))] == FMMHeap<T>::FMM_FAR )	{
				sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i-1,j);
				lv.phi(i-1,j) = nearPhi;
				heap->push(nearS, nearPhi);
			}
			if( i != lv.phi.nX-1 && sTag2D[(nearS = lv.iS(i+1,j))] == FMMHeap<T>::FMM_FAR )	{
				sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i+1,j);
				lv.phi(i+1,j) = nearPhi;
				heap->push(nearS, nearPhi);
			}
			if( j != 0 && sTag2D[(nearS = lv.iS(i,j-1))] == FMMHeap<T>::FMM_FAR )	{
				sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j-1);
				lv.phi(i,j-1) = nearPhi;
				heap->push(nearS, nearPhi);
			}
			if( j != lv.phi.nY-1 && sTag2D[(nearS = lv.iS(i,j+1))] == FMMHeap<T>::FMM_FAR )	{
				sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j+1);
				lv.phi(i,j+1) = nearPhi;
				heap->push(nearS, nearPhi);
			}
		}

//		lv.phi.negate();
//		continue;

		while(!heap->empty())	{
			NodeHeap<T> nodeHeap = heap->pop();

			grid_index_type i = nodeHeap.idx / sTag2D.nY;
			grid_index_type j = nodeHeap.idx % sTag2D.nY;

			sTag2D[nodeHeap.idx] = FMMHeap<T>::FMM_ACCEPTED;
			lv.phi[nodeHeap.idx] = nodeHeap.distance;
			if(nodeHeap.distance > narrowBandWidth)
				break;

			grid_index_type nearS;
			T nearPhi;

			if( i != 0 )	{
				nearS = lv.iS(i-1,j);
				if(sTag2D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0f;
					heap->push(nearS, nearPhi);
				}
				else if(sTag2D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i-1, j);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( i != lv.phi.nX-1 )	{
				nearS = lv.iS(i+1,j);
				if(sTag2D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0f;
					heap->push(nearS, nearPhi);
				}
				else if(sTag2D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i+1, j);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( j != 0 )	{
				nearS = lv.iS(i,j-1);
				if(sTag2D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0f;
					heap->push(nearS, nearPhi);
				}
				else if(sTag2D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j-1);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( j != lv.phi.nY-1 )	{
				nearS = lv.iS(i,j+1);
				if(sTag2D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag2D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0f;
					heap->push(nearS, nearPhi);
				}
				else if(sTag2D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j+1);
					heap->assignDistance(nearS, nearPhi);
				}
			}
		}

		lv.phi.negate();
	}

	lv.phi.scale(lv.h);
	lv.phi.setNeumannBoundary();
}

template<typename T>
void Levelset_RedistancingFMM<T>::redistance_initialize(Levelset2D<T>& lv) const {
	
	interfaceCellList.clear();

	// Data Structure Initialization
	grid_index_type _narrowBandWidth = MIN(static_cast<T>(MAX(lv.nX, lv.nY)), narrowBandWidth);
	minHeapSize = static_cast<grid_index_type>(3.0f * _narrowBandWidth * MAX(lv.nX,lv.nY));

	if(heap == NULL)	heap = new FMMHeap<T>(minHeapSize, lv.numS);
	else				heap->init(minHeapSize);

	if(sTag2D.size() == 0) sTag2D.resize(lv.phi);
	sTag2D.set(FMMHeap<T>::FMM_FAR);

	Grid2<T,T> phi0(lv.phi);

	T extreme = narrowBandWidth;
	for(grid_index_type i=0;i<lv.phi.nX;++i)	{
		for(grid_index_type j=0;j<lv.phi.nY;++j)	{

			int sign0 = SIGN(lv.phi(i,j));

			if( (i != (lv.phi.nX-1)		&& sign0 != SIGN(lv.phi(i+1,j))) 
			 || (j != (lv.phi.nY-1)		&& sign0 != SIGN(lv.phi(i,j+1))) 
			 || (i != 0					&& sign0 != SIGN(lv.phi(i-1,j))) 
			 || (j != 0					&& sign0 != SIGN(lv.phi(i,j-1))) )	{
					interfaceCellList.push_back(i);
					interfaceCellList.push_back(j);
			}
			else
				lv.phi(i,j) = sign0 * extreme;
		}
	}
	//cout << interfaceCellList.size() << endl;

#define POINTRATIO(a, b)	(abs(a/(a-b)))

	// Distance, sTag2D Initialization near interface
	T s, t;

	for(size_t l=0;l<interfaceCellList.size();l+=2)	{
		grid_index_type i = interfaceCellList[l];
		grid_index_type j = interfaceCellList[l+1];
		grid_index_type cell = lv.iS(i,j);

		int sign0 = SIGN(lv.phi[cell]);

		s = t = FMMHeap<T>::MAX_DISTANCE;

		if(i != (lv.phi.nX-1)	&& SIGN(phi0[lv.iS(i+1,j)]) != sign0) s = MIN(s, POINTRATIO(phi0[cell], phi0[lv.iS(i+1,j)]));
		if(i != 0				&& SIGN(phi0[lv.iS(i-1,j)]) != sign0) s = MIN(s, POINTRATIO(phi0[cell], phi0[lv.iS(i-1,j)]));
		if(j != (lv.phi.nY-1)	&& SIGN(phi0[lv.iS(i,j+1)]) != sign0) t = MIN(t, POINTRATIO(phi0[cell], phi0[lv.iS(i,j+1)]));
		if(j != 0				&& SIGN(phi0[lv.iS(i,j-1)]) != sign0) t = MIN(t, POINTRATIO(phi0[cell], phi0[lv.iS(i,j-1)]));

		if(s == FMMHeap<T>::MAX_DISTANCE)		lv.phi[cell] = sign0 * t;
		else if(t == FMMHeap<T>::MAX_DISTANCE)	lv.phi[cell] = sign0 * s;
		else {
			T denominator = sqrt(s*s + t*t);

			if(IsAlmostZero(denominator))	lv.phi[cell] = 0.0f;
			else							lv.phi[cell] = sign0 * s * t / denominator;
		}

		sTag2D[cell] = FMMHeap<T>::FMM_ACCEPTED;

	}

#undef POINTRATIO
}

template<typename T>
const T Levelset_RedistancingFMM<T>::recomputeDistance(const Levelset2D<T>& lv, grid_index_type i, grid_index_type j) const {
	T s, t;
	s = t = FMMHeap<T>::MAX_DISTANCE;

	if( i > 0		&& sTag2D[lv.iS(i-1,j)] == FMMHeap<T>::FMM_ACCEPTED ) s = MIN(s, lv.phi[lv.iS(i-1,j)]);
	if( i < lv.nX+1 && sTag2D[lv.iS(i+1,j)] == FMMHeap<T>::FMM_ACCEPTED ) s = MIN(s, lv.phi[lv.iS(i+1,j)]);
	if( j > 0		&& sTag2D[lv.iS(i,j-1)] == FMMHeap<T>::FMM_ACCEPTED ) t = MIN(t, lv.phi[lv.iS(i,j-1)]);
	if( j < lv.nY+1 && sTag2D[lv.iS(i,j+1)] == FMMHeap<T>::FMM_ACCEPTED ) t = MIN(t, lv.phi[lv.iS(i,j+1)]);

	if( s >= t + 1 )		return t + 1.0f;
	else if( t >= s + 1)	return s + 1.0f;
	else					return 0.5f * ( s + t + sqrt(2.0f - SQ(s-t)) );
}

};		// fluidspace