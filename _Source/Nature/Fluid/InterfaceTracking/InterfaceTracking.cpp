#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking.h>

using namespace std;
using namespace fluidspace;

InterfaceTracking::InterfaceTracking() {
	eye = Vector3f(0, 0, 2);
	ref = Vector3f(0, 0, 0);	
	up  = Vector3f(0.0f, 1.0f, 0.0f);
	theta = 45.0f;	
}