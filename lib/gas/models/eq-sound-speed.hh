// Author: Daniel F Potter
// Version:
//   21-Sep-2009
//      Initial coding.
//   23-Sep-2009
//      Some refactoring to keep coding style consistent
//      with rest of module. (RJG)
//

#ifndef EQ_SOUND_SPEED_HH
#define EQ_SOUND_SPEED_HH

#include "gas_data.hh"
#include "gas-model.hh"
#include "sound-speed-model.hh"

class Eq_sound_speed_model : public Sound_speed_model {
public:
    Eq_sound_speed_model(Gas_model &gm);
    ~Eq_sound_speed_model();

private:
    Gas_model &gm_;
    int s_eval_sound_speed(Gas_data &Q);
};

#endif
