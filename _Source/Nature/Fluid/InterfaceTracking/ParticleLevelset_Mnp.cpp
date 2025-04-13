#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/ParticleLevelset_Mnp.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

PLS_ParticleMnp::PLS_ParticleMnp() {
	r_min_frac = b_min_frac = 0.1f;
	r_max_frac = 0.5f;
	b_max_frac = 3.0f;

	numParticlesPerCell_1D = 4;
	reseed_cycle = 20;
}

void PLS_ParticleMnp::initialize(InterfaceTracking_PLS2D& pls) const {
	
	real h = pls.h;

	r_min = r_min_frac * h;	
	b_min = b_min_frac * h;
	r_max = r_max_frac * h;	
	b_max = b_max_frac * h;
	numParticlesPerCell = numParticlesPerCell_1D * numParticlesPerCell_1D;
	
	pls.particles.clear();
	pls.idxEscapedParticles.clear();

	//////////////////////////////////////////////////////////////////////////
	// For Attracting Particles
	Grid2<Vector2<real>, real> normals;
	normals.resize(pls.phi);
	MacGridLibrary<real>::getNormals(normals, pls.phi);

	//////////////////////////////////////////////////////////////////////////
	// Initial Particle Seeding
	for(grid_index_type i=1;i<pls.phi.nX-1;++i) 
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) 
			add_particles(pls, Vector2<grid_index_type>(i,j), numParticlesPerCell, normals);
}

void PLS_ParticleMnp::initialize(InterfaceTracking_PLS3D& pls) const {

	real h = pls.h;

	r_min = r_min_frac * h;	
	b_min = b_min_frac * h;
	r_max = r_max_frac * h;	
	b_max = b_max_frac * h;
	numParticlesPerCell = numParticlesPerCell_1D * numParticlesPerCell_1D * numParticlesPerCell_1D;

	pls.particles.clear();
	pls.idxEscapedParticles.clear();

	//////////////////////////////////////////////////////////////////////////
	// For Attracting Particles
	Grid3<Vector3<real>, real> normals;
	normals.resize(pls.phi);
	MacGridLibrary<real>::getNormals(normals, pls.phi);

	//////////////////////////////////////////////////////////////////////////
	// Initial Particle Seeding
	for(grid_index_type i=1;i<pls.phi.nX-1;++i) 
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) 
			for(grid_index_type k=1;k<pls.phi.nZ-1;++k) 
				add_particles(pls, Vector3<grid_index_type>(i,j,k), numParticlesPerCell, normals);
}

void PLS_ParticleMnp::reseed(InterfaceTracking_PLS2D& pls) const {
	typedef std::pair<InterfaceTracking_PLS2D::Particle *, real> particlePriority;

	struct compareParticlePriority {
		const bool operator()(const particlePriority& p0, const particlePriority& p1) const {
			return (p0.second < p1.second);
		} 
	};

	typedef std::priority_queue<particlePriority, vector<particlePriority>, compareParticlePriority> particleHeap;

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	log_info << "Before Reseeding" << endl;
	pls.printParticleInformation(log_info);

	real h = pls.h;

	Grid2<particleHeap, real> heaps;
	heaps.resize(pls.phi);

	collectParticles<InterfaceTracking_PLS2D, Grid2<particleHeap,real>, Vector2<grid_index_type>>(heaps, pls);

	//////////////////////////////////////////////////////////////////////////
	// For Attracting Particles
	Grid2<Vector2<real>, real> normals;
	normals.resize(pls.phi);
	MacGridLibrary<real>::getNormals(normals, pls.phi);

	//////////////////////////////////////////////////////////////////////////
	// Reseed Particles for unresolved region
	for(grid_index_type i=1;i<pls.phi.nX-1;++i) 
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) 
			add_particles(pls, Vector2<grid_index_type>(i,j), numParticlesPerCell - heaps(i,j).size(), normals);

	//////////////////////////////////////////////////////////////////////////
	// Delete deleted particles
	delete_isDeletedParticles(pls);

	log_info << "After Reseeding" << endl;
	pls.printParticleInformation(log_info);
}

void PLS_ParticleMnp::reseed(InterfaceTracking_PLS3D& pls) const {

	typedef std::pair<InterfaceTracking_PLS3D::Particle *, real> particlePriority;

	struct compareParticlePriority {
		const bool operator()(const particlePriority& p0, const particlePriority& p1) const {
			return (p0.second < p1.second);
		} 
	};

	typedef std::priority_queue<particlePriority, vector<particlePriority>, compareParticlePriority> particleHeap;

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	log_info << "Before Reseeding" << endl;
	pls.printParticleInformation(log_info);

	real h = pls.h;

	Grid3<particleHeap, real> heaps;
	heaps.resize(pls.phi);

	collectParticles<InterfaceTracking_PLS3D, Grid3<particleHeap,real>, Vector3<grid_index_type>>(heaps, pls);

	log_info << "Collect Particle Done" << endl;

	//////////////////////////////////////////////////////////////////////////
	// For Attracting Particles
	Grid3<Vector3<real>, real> normals;
	normals.resize(pls.phi);
	MacGridLibrary<real>::getNormals(normals, pls.phi);

	log_info << "get Normal Done" << endl;

	//////////////////////////////////////////////////////////////////////////
	// Reseed Particles for unresolved region
	for(grid_index_type i=1;i<pls.phi.nX-1;++i) 
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) 
			for(grid_index_type k=1;k<pls.phi.nZ-1;++k) 
				add_particles(pls, Vector3<grid_index_type>(i,j,k), numParticlesPerCell - heaps(i,j,k).size(), normals);

	log_info << "Add Particles Done" << endl;

	//////////////////////////////////////////////////////////////////////////
	// Delete deleted particles
	delete_isDeletedParticles(pls);

	log_info << "After Reseeding" << endl;
	pls.printParticleInformation(log_info);
}
