// Author: Rowan J. Gollan
// Date: 15-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#ifndef THERMO_INTERPOLATOR_HH
#define THERMO_INTERPOLATOR_HH

#include <string>

#include "../../../lib/gas/models/gas_data.hh"

class Thermo_interpolator {
public:
    Thermo_interpolator();
    virtual ~Thermo_interpolator();

    int one_d_interp(Gas_data &gL1, Gas_data &gL0,
		     Gas_data &gR0, Gas_data &gR1,
		     double cL1Length, double cL0Length,
		     double cR0Length, double cR1Length,
		     Gas_data &Lft, Gas_data &Rght)
    { return s_one_d_interp(gL1, gL0, gR0, gR1,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft, Rght); }
private:
    virtual int s_one_d_interp(Gas_data &gL1, Gas_data &gL0,
			       Gas_data &gR0, Gas_data &gR1,
			       double cL1Length, double cL0Length,
			       double cR0Length, double cR1Length,
			       Gas_data &Lft, Gas_data &Rght) = 0;

};

Thermo_interpolator* create_Thermo_interpolator(std::string name);

#endif
