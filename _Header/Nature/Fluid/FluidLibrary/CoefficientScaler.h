#pragma once

#include "../../Nature.h"
#include "../Environment/Environments.h"

namespace fluidspace {

template<typename T = fluidspace::real>
class Coefficient_Scaler {
public : 
	Coefficient_Scaler(const T original_dx, const T target_dx) : src(original_dx), des(target_dx), factor(des / src) {}

	const T gravity(const T g) const						{ return g * factor;						}
	const T tension_sigma(const T sigma) const				{ return sigma;								}
	const T tension_sigma_freeSurface(const T sigma) const	{ return sigma * (factor*factor*factor);	}

	const T density(const T den) const						{ return den / (factor*factor*factor);		}
	const T viscosity(const T visc) const					{ return visc / factor;						}

private : 
	const T src, des, factor;
};

};		// fluidspace