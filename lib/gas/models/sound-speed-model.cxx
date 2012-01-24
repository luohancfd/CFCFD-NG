// Author: Daniel F Potter
// Version: 
//   21-Sep-2009
//      Initial coding.
//   23-Sep-2009
//       Some refactoring to keep style consistent
//       with rest of module. (RJG)

#include <string>

#include "../../util/source/lua_service.hh"
#include "sound-speed-model.hh"
#include "eq-sound-speed.hh"
#include "noneq-sound-speed.hh"

using namespace std;

Sound_speed_model* create_sound_speed_model(string model, Gas_model &gm, lua_State *L)
{
    if ( model == "equilibrium" ) {
	return new Eq_sound_speed_model(gm);
    }
    else if ( model == "nonequilibrium" ) {
	return new Noneq_sound_speed_model(gm, L);
    }
    else {
	ostringstream ost;
	ost << "create_sound_speed_model():\n";
	ost << "The 'sound_speed': " << model << endl;
	ost << "is unknown or not implemented.\n";
	input_error(ost);
    }
    // If we get here, failed to create a model
    return 0;

}
