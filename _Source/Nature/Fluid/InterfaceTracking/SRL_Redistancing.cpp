#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/SRL_Redistancing.h>

using namespace std;
using namespace fluidspace;

const bool SRL_Redistancing::is_times_up(real dt, real C_CFL) const { 

	cfl_acc += (dt * C_CFL);
	num_updates += 1;

	if(num_updates < min_updates) 
		return false;

	if(num_updates >= max_updates) {
		log_info << "SRL_Redistancing is Needed _ Max Updates Satisfied _ " << max_updates << " Updates Done _ CFL is " << cfl_acc << endl;

		cfl_acc = 0;
		num_updates = 0;		

		return true;
	}

	if(cfl_acc > cfl_frequency) {
		log_info << "SRL_Redistancing is Needed _ CFL Update Criteria is Satisfied _ " << cfl_acc << " cfl _ Number of updates is " << num_updates << endl;

		cfl_acc = 0;
		num_updates = 0;	

		return true;
	}

	return false;
}

void SRL_Redistancing::importBinary(std::ifstream& ifs) {
	::importBinary(ifs, min_updates);
	::importBinary(ifs, max_updates);
	::importBinary(ifs, cfl_frequency);
	::importBinary(ifs, cfl_acc);
	::importBinary(ifs, num_updates);
}

void SRL_Redistancing::exportBinary(std::ofstream& ofs) const {
	::exportBinary(ofs, min_updates);
	::exportBinary(ofs, max_updates);
	::exportBinary(ofs, cfl_frequency);
	::exportBinary(ofs, cfl_acc);
	::exportBinary(ofs, num_updates);
}

void SRL_Redistancing_PDE::post_process() const {
	log_info << std::endl << "### Log for SRL_Redistancing_PDE _ Post Process ###" << std::endl;

	log_info << "Time for total spent time redistancing process             _ " << time_total << std::endl;
	log_info << "Time for total spent time only for fine grid redistancing  _ " << time_only_fine_grid << std::endl << std::endl;
}