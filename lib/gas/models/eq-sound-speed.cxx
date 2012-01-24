// Author: Daniel F Potter
// Version:
//   21-Sep-2009
//      Initial coding.
//   23-Sep-2009
//      Some refactoring to keep coding style consistent
//      with rest of module. (RJG)
//

#include <cmath>

#include "eq-sound-speed.hh"

Eq_sound_speed_model::
Eq_sound_speed_model(Gas_model &gm)
    : gm_(gm) {}

Eq_sound_speed_model::
~Eq_sound_speed_model() {}

int
Eq_sound_speed_model::
s_eval_sound_speed(Gas_data &Q)
{
    // Reference:
    // Cengel and Boles (1998)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // McGraw Hill
    // Equation 16-10 on p. 849
    
    // "frozen" sound speed
    
    int status;
    Q.a = sqrt(gm_.gamma(Q, status)*gm_.dpdrho_const_T(Q, status));
    return status;
}
