// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef SOUND_SPEED_MODEL_HH
#define SOUND_SPEED_MODEL_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas-model.hh"
#include "gas_data.hh"

class Sound_speed_model {
public:
    Sound_speed_model() {};
    virtual ~Sound_speed_model() {}

    int eval_sound_speed(Gas_data &Q)
    { return s_eval_sound_speed(Q); }
    
private:
    virtual int s_eval_sound_speed(Gas_data &Q) = 0;
};

Sound_speed_model* create_sound_speed_model(std::string model, Gas_model &gm, lua_State *L);

#endif
