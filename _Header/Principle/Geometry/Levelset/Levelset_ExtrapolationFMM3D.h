#pragma once

namespace geometryspace {

template<typename T>
void Levelset_ExtrapolationFMM<T>::extrapolate(Grid3<T,T>& val, const GridBoundaryCondition<T,T>& bdry, const Grid3<T,T>& phi) const {
	basis::throwError("Not Implemented _ Levelset_ExtrapolationFMM<T>::extrapolate()");
}

template<typename T>
void Levelset_ExtrapolationFMM<T>::extrapolate(Grid3<T,T> &val, typename MacGrid3D<T>::FIELD_TYPE _type, const GridBoundaryCondition<T,T>& bdry, const Grid3<T,T>& phi) const	{
	
	grid_index_type _narrowBandWidth = MIN(sqrt((real)(phi.nX*phi.nX+phi.nY*phi.nY+phi.nZ*phi.nZ)), narrowBandWidth);
	minHeapSize = static_cast<grid_index_type>(3.0f * narrowBandWidth * SQ(MAX3(phi.nX,phi.nY,phi.nZ)));

	if(heap == NULL)	heap = new FMMHeap<T>(minHeapSize, phi.size());
	else				heap->init(minHeapSize);

	grid_index_type idx[3][2];
	Grid3<char,T> valTag(val.nX, val.nY, val.nZ, val.dx, val.dy, val.dz, val.offset);

	switch(_type)	{
		case MacGrid3D<T>::TYPE_U :	idx[0][0] = 0; idx[0][1] = 1;
									idx[1][0] = 0; idx[1][1] = 0;
									idx[2][0] = 0; idx[2][1] = 0;
									break;
		case MacGrid3D<T>::TYPE_V :	idx[0][0] = 0; idx[0][1] = 0;
									idx[1][0] = 0; idx[1][1] = 1;
									idx[2][0] = 0; idx[2][1] = 0;
									break;
		case MacGrid3D<T>::TYPE_W :	idx[0][0] = 0; idx[0][1] = 0;
									idx[1][0] = 0; idx[1][1] = 0;
									idx[2][0] = 0; idx[2][1] = 1;
									break;
	}

#define SofVEL(i,j,k) (0.5*(phi(i+idx[0][0],j+idx[1][0],k+idx[2][0])+phi(i+idx[0][1],j+idx[1][1],k+idx[2][1])))

	for(grid_index_type i=0;i<val.nX;++i)	{
		for(grid_index_type j=0;j<val.nY;++j)	{
			for(grid_index_type k=0;k<val.nZ;++k)	{
				if(SofVEL(i,j,k) < 0.0f)	valTag(i,j,k) = FMMHeap<T>::FMM_ACCEPTED;
				else						valTag(i,j,k) = FMMHeap<T>::FMM_FAR;
			}			
		}
	}

	for(grid_index_type i=0;i<val.nX;++i)	{
		for(grid_index_type j=0;j<val.nY;++j)	{
			for(grid_index_type k=0;k<val.nZ;++k)	{
				if(valTag(i,j,k)==FMMHeap<T>::FMM_ACCEPTED)	{
					if(i != 0 && valTag(i-1,j,k) == FMMHeap<T>::FMM_FAR)	{
						valTag(i-1,j,k) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i-1,j,k), SofVEL(i-1,j,k));
					}
					if(i != val.nX-1 && valTag(i+1,j,k) == FMMHeap<T>::FMM_FAR)	{
						valTag(i+1,j,k) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i+1,j,k), SofVEL(i+1,j,k));
					}
					if(j != 0 && valTag(i,j-1,k) == FMMHeap<T>::FMM_FAR)	{
						valTag(i,j-1,k) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i,j-1,k), SofVEL(i,j-1,k));
					}
					if(j != val.nY-1 && valTag(i,j+1,k) == FMMHeap<T>::FMM_FAR)	{
						valTag(i,j+1,k) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i,j+1,k), SofVEL(i,j+1,k));
					}
					if(k != 0 && valTag(i,j,k-1) == FMMHeap<T>::FMM_FAR)	{
						valTag(i,j,k-1) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i,j,k-1), SofVEL(i,j,k-1));
					}
					if(k != val.nZ-1 && valTag(i,j,k+1) == FMMHeap<T>::FMM_FAR)	{
						valTag(i,j,k+1) = FMMHeap<T>::FMM_TRIAL;
						heap->push(val.idx(i,j,k+1), SofVEL(i,j,k+1));
					}
				}
			}			
		}
	}

	while(!heap->empty())	{
		NodeHeap<T> nodeHeap = heap->pop();

		grid_index_type i = nodeHeap.idx / (val.nY*val.nZ);
		grid_index_type j = (nodeHeap.idx - i * val.nY * val.nZ) / val.nZ;
		grid_index_type k = nodeHeap.idx % val.nZ;

		if(SofVEL(i,j,k) > narrowBandWidth * phi.dx)
			break;

		valTag(nodeHeap.idx) = FMMHeap<T>::FMM_ACCEPTED;		

		real s,t,u, val_s, val_t, val_u;
		s = t = u = FMMHeap<T>::MAX_DISTANCE;

		if(i != 0 && valTag(i-1,j,k) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i-1,j,k) < s)	{
				val_s = val(i-1,j,k);
				s = SofVEL(i-1,j,k);
			}
		}
		if(i != val.nX-1 && valTag(i+1,j,k) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i+1,j,k) < s)	{
				val_s = val(i+1,j,k);
				s = SofVEL(i+1,j,k);
			}
		}
		if(j != 0 && valTag(i,j-1,k) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i,j-1,k) < t)	{
				val_t = val(i,j-1,k);
				t = SofVEL(i,j-1,k);
			}
		}
		if(j != val.nY-1 && valTag(i,j+1,k) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i,j+1,k) < t)	{
				val_t = val(i,j+1,k);
				t = SofVEL(i,j+1,k);
			}
		}
		if(k != 0 && valTag(i,j,k-1) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i,j,k-1) < u)	{
				val_u = val(i,j,k-1);
				u = SofVEL(i,j,k-1);
			}
		}
		if(k != val.nZ-1 && valTag(i,j,k+1) == FMMHeap<T>::FMM_ACCEPTED)	{
			if(SofVEL(i,j,k+1) < u)	{
				val_u = val(i,j,k+1);
				u = SofVEL(i,j,k+1);
			}
		}

		if(s==FMMHeap<T>::MAX_DISTANCE && t==FMMHeap<T>::MAX_DISTANCE)		val(nodeHeap.idx) = val_u;
		else if(t==FMMHeap<T>::MAX_DISTANCE && u==FMMHeap<T>::MAX_DISTANCE)	val(nodeHeap.idx) = val_s;
		else if(u==FMMHeap<T>::MAX_DISTANCE && s==FMMHeap<T>::MAX_DISTANCE)	val(nodeHeap.idx) = val_t;
		else if(s==FMMHeap<T>::MAX_DISTANCE)	{
			real wt = abs(SofVEL(i,j,k)-t);
			real wu = abs(SofVEL(i,j,k)-u);

			if(wt+wu < 1.e-5)	val(nodeHeap.idx) = 0.5 * (val_t+val_u);
			else				val(nodeHeap.idx) = (wt*val_t + wu*val_u) / (wt+wu);
		}
		else if(t==FMMHeap<T>::MAX_DISTANCE)	{
			real ws = abs(SofVEL(i,j,k)-s);
			real wu = abs(SofVEL(i,j,k)-u);

			if(ws+wu < 1.e-5)	val(nodeHeap.idx) = 0.5 * (val_s+val_u);
			else				val(nodeHeap.idx) = (ws*val_s + wu*val_u) / (ws+wu);
		}
		else if(u==FMMHeap<T>::MAX_DISTANCE)	{
			real wt = abs(SofVEL(i,j,k)-t);
			real ws = abs(SofVEL(i,j,k)-s);

			if(wt+ws < 1.e-5)	val(nodeHeap.idx) = 0.5 * (val_t+val_s);
			else				val(nodeHeap.idx) = (wt*val_t + ws*val_s) / (wt+ws);
		}
		else {
			real ws = abs(SofVEL(i,j,k)-s);
			real wt = abs(SofVEL(i,j,k)-t);
			real wu = abs(SofVEL(i,j,k)-u);			

			if(ws+wt+wu < 1.e-5)	val(nodeHeap.idx) = (val_s+val_t+val_u) / 3.0f;
			else					val(nodeHeap.idx) = (ws*val_s + wt*val_t + wu*val_u) / (ws+wt+wu);
		}		

		if(i != 0 && valTag(i-1,j,k) == FMMHeap<T>::FMM_FAR)	{
			valTag(i-1,j,k) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i-1,j,k), SofVEL(i-1,j,k));
		}
		if(i != val.nX-1 && valTag(i+1,j,k) == FMMHeap<T>::FMM_FAR)	{
			valTag(i+1,j,k) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i+1,j,k), SofVEL(i+1,j,k));
		}
		if(j != 0 && valTag(i,j-1,k) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j-1,k) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j-1,k), SofVEL(i,j-1,k));
		}
		if(j != val.nY-1 && valTag(i,j+1,k) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j+1,k) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j+1,k), SofVEL(i,j+1,k));
		}
		if(k != 0 && valTag(i,j,k-1) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j,k-1) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j,k-1), SofVEL(i,j,k-1));
		}
		if(k != val.nZ-1 && valTag(i,j,k+1) == FMMHeap<T>::FMM_FAR)	{
			valTag(i,j,k+1) = FMMHeap<T>::FMM_TRIAL;
			heap->push(val.idx(i,j,k+1), SofVEL(i,j,k+1));
		}
	}

	for(grid_index_type i=0;i<val.size();i++)
		if(valTag(i) != FMMHeap<T>::FMM_ACCEPTED)	val(i) = 0.0f;

#undef SofVEL

	bdry.enforce(val);
}

};			// geometryspace