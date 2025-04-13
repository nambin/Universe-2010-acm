#pragma once

namespace geometryspace {

template<typename T>
void Levelset_RedistancingFMM<T>::redistance(Levelset3D<T>& lv)	const {

	// Fast Marching in New NarrowBand
//	lv.phi.setNeumannBoundary();
	lv.phi.scale(lv.hInv);

	// Data Structure & 
	// Distance, cTag Initialization near interface
	redistance_initialize(lv);

	// FMM initialization for positive, negative Region
	for(int k=0;k<2;k++)	{
		if(k==1)	
			heap->init(minHeapSize);

		for(size_t l=0;l<interfaceCellList.size();l+=3)	{
			grid_index_type i = interfaceCellList[l];
			grid_index_type j = interfaceCellList[l+1];
			grid_index_type k = interfaceCellList[l+2];

			if( lv.phi(i,j,k) <= 0 )
				continue;		

			grid_index_type nearS;
			T nearPhi;			

			if( i != 0 && sTag3D[(nearS = lv.iS(i-1,j,k))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i-1,j,k);
				heap->push(nearS, nearPhi);
			}
			if( i != lv.phi.nX-1 && sTag3D[(nearS = lv.iS(i+1,j,k))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i+1,j,k);
				heap->push(nearS, nearPhi);
			}
			if( j != 0 && sTag3D[(nearS = lv.iS(i,j-1,k))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j-1,k);
				heap->push(nearS, nearPhi);
			}
			if( j != lv.phi.nY-1 && sTag3D[(nearS = lv.iS(i,j+1,k))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j+1,k);
				heap->push(nearS, nearPhi);
			}
			if( k != 0 && sTag3D[(nearS = lv.iS(i,j,k-1))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j,k-1);
				heap->push(nearS, nearPhi);
			}
			if( k != lv.phi.nZ-1 && sTag3D[(nearS = lv.iS(i,j,k+1))] == FMMHeap<T>::FMM_FAR )	{
				sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
				nearPhi = recomputeDistance(lv, i,j,k+1);
				heap->push(nearS, nearPhi);
			}
		}

		while(!heap->empty())	{
			NodeHeap<T> nodeHeap = heap->pop();

			grid_index_type i = nodeHeap.idx / (sTag3D.nY*sTag3D.nZ);
			grid_index_type j = (nodeHeap.idx - i * sTag3D.nY * sTag3D.nZ) / sTag3D.nZ;
			grid_index_type k = nodeHeap.idx % sTag3D.nZ;

			sTag3D[nodeHeap.idx] = FMMHeap<T>::FMM_ACCEPTED;
			lv.phi[nodeHeap.idx] = nodeHeap.distance;
			if(nodeHeap.distance > narrowBandWidth)
				break;

			grid_index_type nearS;
			T nearPhi;

			if( i != 0 )	{
				nearS = lv.iS(i-1,j,k);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i-1, j, k);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( i != lv.phi.nX-1 )	{
				nearS = lv.iS(i+1,j,k);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i+1, j, k);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( j != 0 )	{
				nearS = lv.iS(i,j-1,k);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j-1, k);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( j != lv.phi.nY-1 )	{
				nearS = lv.iS(i,j+1,k);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j+1, k);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( k != 0 )	{
				nearS = lv.iS(i,j,k-1);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j, k-1);
					heap->assignDistance(nearS, nearPhi);
				}
			}
			if( k != lv.phi.nZ-1 )	{
				nearS = lv.iS(i,j,k+1);
				if(sTag3D[nearS] == FMMHeap<T>::FMM_FAR)	{
					sTag3D[nearS] = FMMHeap<T>::FMM_TRIAL;
					nearPhi = lv.phi[nodeHeap.idx] + 1.0;
					heap->push(nearS, nearPhi);
				}
				else if(sTag3D[nearS] == FMMHeap<T>::FMM_TRIAL)	{
					nearPhi = recomputeDistance(lv, i, j, k+1);
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
void Levelset_RedistancingFMM<T>::redistance_initialize(Levelset3D<T>& lv) const {

	interfaceCellList.clear();

	T maxDist = MAX3(lv.nX, lv.nY, lv.nZ);

	// Data Structure Initialization
	grid_index_type _narrowBandWidth = MIN(sqrt((T)(lv.nX*lv.nX+lv.nY*lv.nY+lv.nZ*lv.nZ)), narrowBandWidth);
	minHeapSize = static_cast<grid_index_type>(3.0f * narrowBandWidth * SQ(MAX3(lv.nX,lv.nY,lv.nZ)));

	if(heap == NULL)	heap = new FMMHeap<T>(minHeapSize, lv.numS);
	else				heap->init(minHeapSize);

	if(sTag3D.size() == 0) sTag3D.resize(lv.phi);
	sTag3D.set(FMMHeap<T>::FMM_FAR);

	Grid3<T,T> phi0(lv.phi);

	for(grid_index_type i=0;i<lv.phi.nX;++i)	{
		for(grid_index_type j=0;j<lv.phi.nY;++j)	{
			for(grid_index_type k=0;k<lv.phi.nZ;++k)	{

				int sign0 = SIGN(lv.phi(i,j,k));
				lv.phi(i,j,k) = sign0 * maxDist; 

				if( (i != (lv.phi.nX-1) && sign0 != SIGN(lv.phi(i+1,j,k))) 
				|| (j != (lv.phi.nY-1) && sign0 != SIGN(lv.phi(i,j+1,k))) 
				|| (k != (lv.phi.nZ-1) && sign0 != SIGN(lv.phi(i,j,k+1))) 
				|| (i != 0          && sign0 != SIGN(lv.phi(i-1,j,k))) 
				|| (j != 0          && sign0 != SIGN(lv.phi(i,j-1,k)))
				|| (k != 0          && sign0 != SIGN(lv.phi(i,j,k-1))) )	{
						interfaceCellList.push_back(i);
						interfaceCellList.push_back(j);
						interfaceCellList.push_back(k);
				}

			}			
		}
	}

#define POINTRATIO(a, b)	(abs(a/(a-b)))

	// Distance, sTag3D Initialization near interface
	T s, t, u;

	for(int l=0;l<(int)interfaceCellList.size();l+=3)	{

		grid_index_type i = interfaceCellList[l];
		grid_index_type j = interfaceCellList[l+1];
		grid_index_type k = interfaceCellList[l+2];
		grid_index_type cell = lv.iS(i,j,k);

		int sign0 = SIGN(lv.phi[cell]);

		s = t = u = FMMHeap<T>::MAX_DISTANCE;

		if(i != (lv.phi.nX-1)	&& SIGN(phi0[lv.iS(i+1,j,k)]) != sign0) s = MIN(s, POINTRATIO(phi0[cell], phi0[lv.iS(i+1,j,k)]));
		if(i != 0				&& SIGN(phi0[lv.iS(i-1,j,k)]) != sign0) s = MIN(s, POINTRATIO(phi0[cell], phi0[lv.iS(i-1,j,k)]));
		if(j != (lv.phi.nY-1)	&& SIGN(phi0[lv.iS(i,j+1,k)]) != sign0) t = MIN(t, POINTRATIO(phi0[cell], phi0[lv.iS(i,j+1,k)]));
		if(j != 0				&& SIGN(phi0[lv.iS(i,j-1,k)]) != sign0) t = MIN(t, POINTRATIO(phi0[cell], phi0[lv.iS(i,j-1,k)]));
		if(k != (lv.phi.nZ-1)	&& SIGN(phi0[lv.iS(i,j,k+1)]) != sign0) u = MIN(u, POINTRATIO(phi0[cell], phi0[lv.iS(i,j,k+1)]));
		if(k != 0				&& SIGN(phi0[lv.iS(i,j,k-1)]) != sign0) u = MIN(u, POINTRATIO(phi0[cell], phi0[lv.iS(i,j,k-1)]));

		if(s == FMMHeap<T>::MAX_DISTANCE && t == FMMHeap<T>::MAX_DISTANCE)		lv.phi[cell] = sign0 * u;
		else if(t == FMMHeap<T>::MAX_DISTANCE && u == FMMHeap<T>::MAX_DISTANCE)	lv.phi[cell] = sign0 * s;
		else if(u == FMMHeap<T>::MAX_DISTANCE && s == FMMHeap<T>::MAX_DISTANCE)	lv.phi[cell] = sign0 * t;
		else if(s == FMMHeap<T>::MAX_DISTANCE)	{
			T denominator = sqrt(SQ(t) + SQ(u));

			if(IsAlmostZero(denominator))	lv.phi[cell] = 0.0f;
			else							lv.phi[cell] = sign0 * t * u / denominator;
		}
		else if(t == FMMHeap<T>::MAX_DISTANCE)	{
			T denominator = sqrt(SQ(u) + SQ(s));

			if(IsAlmostZero(denominator))	lv.phi[cell] = 0.0f;
			else							lv.phi[cell] = sign0 * u * s / denominator;
		}
		else if(u == FMMHeap<T>::MAX_DISTANCE)	{
			T denominator = sqrt(SQ(s) + SQ(t));

			if(IsAlmostZero(denominator))	lv.phi[cell] = 0.0f;
			else							lv.phi[cell] = sign0 * s * t / denominator;
		}
		else	{
			T denominator = sqrt(SQ(s*t) + SQ(t*u) + SQ(u*s));

			if(IsAlmostZero(denominator))	lv.phi[cell] = 0.0f;
			else							lv.phi[cell] = sign0 * s * t * u / denominator;
		}							

		sTag3D[cell] = FMMHeap<T>::FMM_ACCEPTED;

	}

#undef POINTRATIO
}

template<typename T>
const T Levelset_RedistancingFMM<T>::recomputeDistance(Levelset3D<T>& lv, grid_index_type i, grid_index_type j, grid_index_type k) const {

	grid_index_type _iS = lv.iS(i,j,k);
	grid_index_type nX = lv.nX, nY = lv.nY, nZ = lv.nZ;

	T s, t, u;

	s = t = u = FMMHeap<T>::MAX_DISTANCE;

	if( i > 0    && sTag3D[lv.iS(i-1,j,k)] == FMMHeap<T>::FMM_ACCEPTED )	s = MIN(s, lv.phi[lv.iS(i-1,j,k)]);
	if( i < nX+1 && sTag3D[lv.iS(i+1,j,k)] == FMMHeap<T>::FMM_ACCEPTED )	s = MIN(s, lv.phi[lv.iS(i+1,j,k)]);
	if( j > 0    && sTag3D[lv.iS(i,j-1,k)] == FMMHeap<T>::FMM_ACCEPTED )	t = MIN(t, lv.phi[lv.iS(i,j-1,k)]);
	if( j < nY+1 && sTag3D[lv.iS(i,j+1,k)] == FMMHeap<T>::FMM_ACCEPTED )	t = MIN(t, lv.phi[lv.iS(i,j+1,k)]);
	if( k > 0    && sTag3D[lv.iS(i,j,k-1)] == FMMHeap<T>::FMM_ACCEPTED )	u = MIN(u, lv.phi[lv.iS(i,j,k-1)]);
	if( k < nZ+1 && sTag3D[lv.iS(i,j,k+1)] == FMMHeap<T>::FMM_ACCEPTED )	u = MIN(u, lv.phi[lv.iS(i,j,k+1)]);

	if( s == FMMHeap<T>::MAX_DISTANCE && t == FMMHeap<T>::MAX_DISTANCE )			lv.phi[_iS] = u + 1.0f;
	else if( t == FMMHeap<T>::MAX_DISTANCE && u == FMMHeap<T>::MAX_DISTANCE )		lv.phi[_iS] = s + 1.0f;
	else if( s == FMMHeap<T>::MAX_DISTANCE && u == FMMHeap<T>::MAX_DISTANCE )		lv.phi[_iS] = t + 1.0f;
	else if( s == FMMHeap<T>::MAX_DISTANCE )	{
		T a = 2.0f;
		T b = u + t;
		T c = SQ(u) + SQ(t) - 1.0f;

		T det = b*b - a*c;
		T val;

		if(det < 0)	val = ( b ) / a;
		else		val = ( b + sqrt(det) ) / a;

		if(lv.phi[_iS] > val)
			lv.phi[_iS] = val;
	}
	else if( t == FMMHeap<T>::MAX_DISTANCE )	{
		T a = 2.0f;
		T b = u + s;
		T c = SQ(u) + SQ(s) - 1.0f;

		T det = b*b - a*c;
		T val;

		if(det < 0)	val = ( b ) / a;
		else		val = ( b + sqrt(det) ) / a;

		if(lv.phi[_iS] > val)	
			lv.phi[_iS] = val;
	}
	else if( u == FMMHeap<T>::MAX_DISTANCE )	{
		T a = 2.0f;
		T b = t + s;
		T c = SQ(t) + SQ(s) - 1.0f;

		T det = b*b - a*c;
		T val;

		if(det < 0)	val = ( b ) / a;
		else		val = ( b + sqrt(det) ) / a;

		if(lv.phi[_iS] > val)	
			lv.phi[_iS] = val;
	}
	else	{
		T a = 3.0f;
		T b = s + t + u;
		T c = SQ(s) + SQ(t) + SQ(u) - 1.0f;

		T det = b*b - a*c;
		T val;

		if(det < 0)	val = ( b ) / a;
		else		val = ( b + sqrt(det) ) / a;

		if(lv.phi[_iS] > val)	
			lv.phi[_iS] = val;
	}

	//cout << lv.phi[_iS] << " ";
	return lv.phi[_iS];
}

};		// fluidspace