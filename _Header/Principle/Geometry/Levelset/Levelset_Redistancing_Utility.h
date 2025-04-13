#pragma once

namespace geometryspace {

template<typename T>
void Levelset_Redistancing<T>::setPhi_asExtreme(Grid2<T,T>& phi) const {

	T h = MAX(phi.dx, phi.dy);
	T extreme = narrowBandWidth * h;

	for(grid_index_type i=0;i<phi.nX;++i)	{
		for(grid_index_type j=0;j<phi.nY;++j)	{

			int sign0 = SIGN(phi(i,j));

			if(!((i != (phi.nX-1)		&& sign0 != SIGN(phi(i+1,j))) 
			  || (j != (phi.nY-1)		&& sign0 != SIGN(phi(i,j+1))) 
			  || (i != 0				&& sign0 != SIGN(phi(i-1,j))) 
			  || (j != 0				&& sign0 != SIGN(phi(i,j-1))) ))	{

				  phi(i,j) = sign0 * extreme;	

			}
		}
	}
}

template<typename T>
void Levelset_Redistancing<T>::setPhi_asExtreme(Grid3<T,T>& phi) const {
}

};		// geometryspace