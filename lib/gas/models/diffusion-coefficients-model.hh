// Author: Rowan J. Gollan
// Date: 09-Jul-2008

#ifndef DIFFUSION_COEFFICIENTS_MODEL_HH
#define DIFFUSION_COEFFICIENTS_MODEL_HH

class Diffusion_coefficients_model {
public:
    Diffusion_coefficients_model() {}
    virtual ~Diffusion_coefficients_model() {}

    int eval_diffusion_coefficients(Gas_data &Q)
    { return s_eval_diffusion_coefficients(Q); }
private:
    virtual int s_eval_diffusion_coefficients(Gas_data &Q) = 0;
};

#endif
