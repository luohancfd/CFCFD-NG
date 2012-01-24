// Author: Daniel F Potter
// Version: 
//   21-Sep-2009
//      Initial coding.
//   23-Sep-2009
//       Some refactoring to keep style consistent
//       with rest of module. (RJG)

#ifndef NONEQ_SOUND_SPEED_HH
#define NONEQ_SOUND_SPEED_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "gas-model.hh"
#include "sound-speed-model.hh"

class Noneq_sound_speed_model : public Sound_speed_model {
public:
    Noneq_sound_speed_model(Gas_model &gm, lua_State *L);
    ~Noneq_sound_speed_model();

private:
    Gas_model &gm_;
    int iT_, iTe_;
    int s_eval_sound_speed(Gas_data &Q);
};

#endif
