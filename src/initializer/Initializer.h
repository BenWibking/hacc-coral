#ifndef Initializer_Header_Included
#define Initializer_Header_Included

#include "TypesAndDefs.h"
#include "Basedata.h"

void init_particles(::real* pos_x, ::real* pos_y, ::real* pos_z, 
		    ::real* vel_x, ::real* vel_y, ::real* vel_z, 
		    Basedata& bdata, const char *tfName, int useWN = 0);

#endif
