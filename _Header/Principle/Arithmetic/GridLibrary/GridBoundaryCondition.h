#ifndef _FLUIDSPACE_CFD_BOUNDARY_CONDITION_H_
#define _FLUIDSPACE_CFD_BOUNDARY_CONDITION_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "Grid2_Functors.h"
#include "Grid3_Functors.h"
#include "GridLibrary.h"

namespace arithmetic {

using namespace basis;

template<typename T, typename R>
class GridBoundaryCondition : public basis::ParallelComputable {
public : 
	GridBoundaryCondition() {}
	virtual void enforce(Grid2<T,R>&) const = 0;
	virtual void enforce(Grid3<T,R>&) const = 0;

	virtual void print_log(Log& log_) const {}
};

template<typename T, typename R>
class GridBoundaryCondition_Neumann : public GridBoundaryCondition<T,R> {
public : 
	GridBoundaryCondition_Neumann() {}
	void enforce(Grid2<T,R>&) const;
	void enforce(Grid3<T,R>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << std::endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_Neumann"); }
};

template<typename T, typename R>
void GridBoundaryCondition_Neumann<T,R>::enforce(Grid2<T,R>& grid) const {
	Grid2_XNeumann<T,R> xNeumann(grid); compute1D_for(xNeumann, 0, grid.nY);
	Grid2_YNeumann<T,R> yNeumann(grid); compute1D_for(yNeumann, 0, grid.nX);	
}

template<typename T, typename R>
void GridBoundaryCondition_Neumann<T,R>::enforce(Grid3<T,R>& grid) const {
	Grid3_XNeumann<T,R> xNeumann(grid); xNeumann(basis::blocked_range2d<grid_index_type>(0, grid.nY, 0, grid.nZ));
	Grid3_YNeumann<T,R> yNeumann(grid); yNeumann(basis::blocked_range2d<grid_index_type>(0, grid.nX, 0, grid.nZ));
	Grid3_ZNeumann<T,R> zNeumann(grid); zNeumann(basis::blocked_range2d<grid_index_type>(0, grid.nX, 0, grid.nY));
}

};		// arithmetic

#endif