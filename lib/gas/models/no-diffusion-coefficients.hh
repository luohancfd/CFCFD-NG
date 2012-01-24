// Author: Rowan J. Gollan
// Date: 16-Jul-2008

#ifndef NO_DIFFUSION_COEFFICIENTS_HH
#define NO_DIFFUSION_COEFFICIENTS_HH

#include "gas_data.hh"
#include "diffusion-coefficients-model.hh"

class No_diffusion_coefficients : public Diffusion_coefficients_model {
public:
    No_diffusion_coefficients() {}
    ~No_diffusion_coefficients() {}
private:
    int s_eval_diffusion_coefficients(Gas_data &Q);
};

#endif
