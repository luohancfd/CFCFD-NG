// bc_ablating.hh

#include "bc.hh"

class AblatingBC : public BoundaryCondition, public ZeroSystem {
public:
    double Twall;
private:
    std::vector<double> mdot;
    std::string filename;
    std::vector<double> TProfile;
    int ncell_for_profile;
    double mdot_total;
    Gas_model *gmodel;
public:
    AblatingBC( Block *bdp, int which_boundary, double Twall, 
		std::vector<double> &mdot, const std::string filename="fstc_temp.txt" );
    AblatingBC( const AblatingBC &bc );
    virtual ~AblatingBC();
    virtual int apply_inviscid( double t );	// sets ghost cell flow conditions
    virtual int apply_viscous( double t ); 	// sets wall temperature (same as FixedTBC)
private:
    int calculate_ghost_cell_flow_state( FV_Cell *cell1, FV_Interface *wall, FV_Cell *cell0 );
// The following are for the zero system
public:
    int f( const std::valarray<double> &y, std::valarray<double> &G );
    int Jac( const std::valarray<double> &y, Valmatrix &dGdy );
    double eval_du0_drho0i( int isp );
private:
    Gas_data *Q;
    int u0_index;
    int max_iterations;
    double tol;
    std::vector<double> cell_massf;
    double cell_rho;
    double cell_un;
    double cell_mass_flux;
    double cell_momentum_flux;
    Vector3 cell_local_vel;
    
    ZeroFinder * zero_solver;
    double f_jac;
    
    std::valarray<double> y_guess;
    std::valarray<double> y_out;
};
