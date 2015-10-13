#include "spinice.h"

void dipassign(double* dip, int m, bool ori) {

	dip[0] = 1.0;
	dip[1] = 1.0;
	dip[2] = 1.0;

	switch(m)
	{
		case 1: dip[0] *= -1; break;
		case 2: dip[1] *= -1; break;
		case 3: dip[2] *= -1; break;
	}

	//normalization and orientation
	for(int x=0; x<3; x++) {
		dip[x] -= (!ori)*2*dip[x];
		dip[x] /= sqrt(3);
	}
	
	return;
}	
