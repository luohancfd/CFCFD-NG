// bc_sliding_t.hh

#include "bc.hh"

class SlidingTBC : public BoundaryCondition {
public:
    double Twall_i;
    double Twall_f;
    double t_i;
    double t_f;
public:
    SlidingTBC( Block *bdp, int which_boundary, double Twall_i, double Twall_f, double t_i, double t_f );
    SlidingTBC( const SlidingTBC &bc );
    SlidingTBC();
    SlidingTBC & operator=(const SlidingTBC &bc);
    virtual ~SlidingTBC();
    // default apply_inviscid() is just to reflect normal velocity
    virtual int apply_viscous( double t ); // sets wall temperature
};
