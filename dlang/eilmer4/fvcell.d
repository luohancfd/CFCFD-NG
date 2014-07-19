/**
 * fvcell.d
 * Finite-volume cell class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcell;

import std.conv;
import std.string;
import geom;
import gasmodel;
import fvcore;
import flowstate;
import conservedquantities;
import fvvertex;
import fvinterface;
import globalconfig;


class FVCell {
public:
    uint id;  // allows us to work out where, in the block, the cell is
    bool fr_reactions_allowed; // if true, will call chemical_increment (also thermal_increment)
    double dt_chem; // acceptable time step for finite-rate chemistry
    double dt_therm; // acceptable time step for thermal relaxation
    bool in_turbulent_zone; // if true, we will keep the turbulence viscosity
    double base_qdot; // base-level of heat addition to cell, W/m**3
    // Geometry
    Vector3[] pos; // Centre x,y,z-coordinates for time-levels, m,m,m
    double[] volume; // Cell volume for time-levels (per unit depth or radian in 2D), m**3
    double[] areaxy; // (x,y)-plane area for time-levels, m**2
    double iLength; // length in the i-index direction
    double jLength; // length in the j-index direction
    double kLength; // length in the k-index direction
    double L_min;   // minimum length scale for cell
    double distance_to_nearest_wall; // for turbulence model correction.
    double half_cell_width_at_wall;  // ditto
    FVCell cell_at_nearest_wall;   // ditto
    // Connections
    FVInterface[] iface;  // references to defining interfaces of cell
    FVVertex[] vtx;  // references to vertices for quad (2D) and hexahedral (3D) cells
    // Flow
    FlowState fs; // Flow properties
    ConservedQuantities[] U;  // Conserved flow quantities for the update stages.
    ConservedQuantities[] dUdt; // Time derivatives for the update stages.
    ConservedQuantities Q; // source (or production) terms
    // Terms for loose-coupling of radiation.
    double Q_rad_org;
    double f_rad_org;
    double Q_rE_rad; // Rate of energy addition to cell via radiation.
    double Q_rE_rad_save; // Presently, the radiation source term is calculated
                          // at the first update stage.  We need to retain that
                          // value for all of the update stages.
    // Data for computing residuals.
    double rho_at_start_of_step, rE_at_start_of_step;
    // [TODO] implicit variables

    this(in GasModel gm, size_t id_init=0)
    {
	id = id_init;
	pos.length = n_time_levels;
	volume.length = n_time_levels;
	fs = new FlowState(gm, 100.0e3, [300.0,], Vector3(0.0,0.0,0.0));
	foreach(i; 0 .. n_time_levels) {
	    U ~= new ConservedQuantities(gm);
	    dUdt ~= new ConservedQuantities(gm);
	}
	Q = new ConservedQuantities(gm);
    }

    void copy_values_from(in FVCell other, uint type_of_copy)
    {
	switch ( type_of_copy ) {
	case copy_flow_data:
	    fs.copy_values_from(other.fs);
	    Q.copy_values_from(other.Q);
	    foreach(i; 0 .. n_time_levels) {
		U[i].copy_values_from(other.U[i]);
		dUdt[i].copy_values_from(other.dUdt[i]);
	    }
	    break;
	case copy_grid_data:
	    foreach(i; 0 .. n_time_levels) {
		pos[i] = other.pos[i];
		volume[i] = other.volume[i];
		areaxy[i] = other.areaxy[i];
	    }
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    break;
	case copy_cell_lengths_only:
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    break;
	case copy_all_data: 
	default:
	    // [TODO] really need to think about what needs to be copied...
	    id = other.id;
	    foreach(i; 0 .. n_time_levels) {
		pos[i] = other.pos[i];
		volume[i] = other.volume[i];
		areaxy[i] = other.areaxy[i];
	    }
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    fs.copy_values_from(other.fs);
	    Q.copy_values_from(other.Q);
	    foreach(i; 0 .. n_time_levels) {
		U[i].copy_values_from(other.U[i]);
		dUdt[i].copy_values_from(other.dUdt[i]);
	    }
	} // end switch
    }

    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
	pos[to_level] = pos[from_level];
	volume[to_level] = volume[from_level];
	areaxy[to_level] = areaxy[from_level];
	// When working over all cells in a block, the following copies
	// will no doubt do some doubled-up work, but it should be otherwise benign.
	foreach(ref face; iface) {
	    if ( face ) face.copy_grid_level_to_level(from_level, to_level);
	}
	foreach(ref v; vtx) {
	    if ( v ) v.copy_grid_level_to_level(from_level, to_level);
	}
    }

    override string toString()
    {
	char[] repr;
	repr ~= "FVCell(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", pos=" ~ to!string(pos);
	repr ~= ", volume=" ~ to!string(volume);
	repr ~= ", areaxy=" ~ to!string(areaxy);
	repr ~= ", iLength=" ~ to!string(iLength);
	repr ~= ", jLength=" ~ to!string(jLength);
	repr ~= ", kLength=" ~ to!string(kLength);
	repr ~= ", L_min=" ~ to!string(L_min);
	repr ~= ", dt_chem=" ~ to!string(dt_chem);
	repr ~= ", dt_therm=" ~ to!string(dt_therm);
	repr ~= ", in_turbulent_zone=" ~ to!string(in_turbulent_zone);
	repr ~= ", fr_reactions_allowed=" ~ to!string(fr_reactions_allowed);
	repr ~= ", fs=" ~ to!string(fs);
	repr ~= ", U=" ~ to!string(U);
	repr ~= ", dUdt=" ~ to!string(dUdt);
	repr ~= ")";
	return to!string(repr);
    }

    const bool point_is_inside(in Vector3 p, int dimensions, int gtl)
    // Returns true if the point p is inside or on the cell surface.
    {
	if ( dimensions == 2 ) {
	    // In 2 dimensions,
	    // we split the x,y-plane into half-planes and check which side p is on.
	    double xA = vtx[1].pos[gtl].x; double yA = vtx[1].pos[gtl].y;
	    double xB = vtx[1].pos[gtl].x; double yB = vtx[2].pos[gtl].y;
	    double xC = vtx[3].pos[gtl].x; double yC = vtx[3].pos[gtl].y;
	    double xD = vtx[0].pos[gtl].x; double yD = vtx[0].pos[gtl].y;
	    // Now, check to see if the specified point is on the
	    // left of (or on) each bloundary line AB, BC, CD and DA.
	    if ((p.x - xB) * (yA - yB) >= (p.y - yB) * (xA - xB) &&
		(p.x - xC) * (yB - yC) >= (p.y - yC) * (xB - xC) &&
		(p.x - xD) * (yC - yD) >= (p.y - yD) * (xC - xD) &&
		(p.x - xA) * (yD - yA) >= (p.y - yA) * (xD - xA)) {
		return true;
	    } else {
		return false;
	    }
	} else {
	    // In 3 dimensions,
	    // the test consists of dividing the 6 cell faces into triangular facets
	    // with outwardly-facing normals and then computing the volumes of the
	    // tetrahedra formed by these facets and the sample point p.
	    // If any of the tetrahedra volumes are positive
	    // (i.e. p is on the positive side of a facet) and we assume a convex cell,
	    // it means that the point is outside the cell and we may say so
	    // without further testing.

	    // North
	    if ( tetrahedron_volume(vtx[2].pos[gtl], vtx[3].pos[gtl], vtx[7].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[7].pos[gtl], vtx[6].pos[gtl], vtx[2].pos[gtl], p) > 0.0 ) return false;
	    // East
	    if ( tetrahedron_volume(vtx[1].pos[gtl], vtx[2].pos[gtl], vtx[6].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[6].pos[gtl], vtx[5].pos[gtl], vtx[1].pos[gtl], p) > 0.0 ) return false;
	    // South
	    if ( tetrahedron_volume(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[5].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[5].pos[gtl], vtx[4].pos[gtl], vtx[0].pos[gtl], p) > 0.0 ) return false;
	    // West
	    if ( tetrahedron_volume(vtx[3].pos[gtl], vtx[0].pos[gtl], vtx[4].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[4].pos[gtl], vtx[7].pos[gtl], vtx[3].pos[gtl], p) > 0.0 ) return false;
	    // Bottom
	    if ( tetrahedron_volume(vtx[1].pos[gtl], vtx[0].pos[gtl], vtx[3].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[3].pos[gtl], vtx[2].pos[gtl], vtx[1].pos[gtl], p) > 0.0 ) return false;
	    // Top
	    if ( tetrahedron_volume(vtx[4].pos[gtl], vtx[5].pos[gtl], vtx[6].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[6].pos[gtl], vtx[7].pos[gtl], vtx[4].pos[gtl], p) > 0.0 ) return false;
	    // If we arrive here, we haven't determined that the point is outside...
	    return true;
	} // end dimensions != 2
    } // end point_is_inside()

    const void copy_values_to_buffer(ref double[] buf, int type_of_copy, int gtl) {
	throw new Error("[TODO] not yet implemented");
    }

    void copy_values_from_buffer(in double buf, int type_of_copy, int gtl) {
	throw new Error("[TODO] not yet implemented");
    }

    void replace_flow_data_with_average(in FVCell[] src) {
	throw new Error("[TODO] not yet implemented");
    }

    void scan_flow_values_from_string(string bufptr) {
	throw new Error("[TODO] not yet implemented");
    }

    const string write_flow_values_to_string() {
	throw new Error("[TODO] not yet implemented");
    }

    void scan_BGK_from_string(string bufptr) {
	throw new Error("[TODO] not yet implemented");
    }

    const string write_BGK_to_string() {
	throw new Error("[TODO] not yet implemented");
    }

    void encode_conserved(int gtl, int ftl, double omegaz, bool with_k_omega) {
	throw new Error("[TODO] not yet implemented");
    }

    void decode_conserved(int gtl, int ftl, double omegaz, bool with_k_omega) {
	throw new Error("[TODO] not yet implemented");
    }

    bool check_flow_data() {
	throw new Error("[TODO] not yet implemented");
    }


/+ [TODO]
    int time_derivatives(size_t gtl, size_t ftl, size_t dimensions, bool with_k_omega);
    int stage_1_update_for_flow_on_fixed_grid(double dt, bool force_euler, bool with_k_omega);
    int stage_2_update_for_flow_on_fixed_grid(double dt, bool with_k_omega);
    int stage_3_update_for_flow_on_fixed_grid(double dt, bool with_k_omega);
    int stage_1_update_for_flow_on_moving_grid(double dt, bool with_k_omega);
    int stage_2_update_for_flow_on_moving_grid(double dt, bool with_k_omega);
    int chemical_increment(double dt, double T_frozen);
    int thermal_increment(double dt, double T_frozen_energy);
    double signal_frequency(size_t dimensions, bool with_k_omega);
    int turbulence_viscosity_zero();
    int turbulence_viscosity_zero_if_not_in_zone();
    int turbulence_viscosity_limit(double factor);
    int turbulence_viscosity_factor(double factor);
    int turbulence_viscosity_k_omega();
    int update_k_omega_properties(double dt);
    int k_omega_time_derivatives(double *Q_rtke, double *Q_romega, double tke, double omega);
    int clear_source_vector();
    int add_inviscid_source_vector(int gtl, double omegaz=0.0);
    int add_viscous_source_vector(bool with_k_omega);
    double calculate_wall_Reynolds_number(int which_boundary);
    int store_rad_scaling_params(void);
    int rescale_Q_rE_rad(void);
    int reset_Q_rad_to_zero(void);
    double rad_scaling_ratio(void) const;
+/
} // end class FVCell


string[] variable_list_for_cell()
{
    // This function needs to be kept consistent with functions
    // FVCell.write_values_to_string, FVCell.scan_values_from_string
    // (found above) and with the corresponding Python functions
    // write_cell_data and variable_list_for_cell
    // that may be found in app/eilmer3/source/e3_flow.py.
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z"];
    list ~= ["rho", "vel.x", "vel.y", "vel.z"];
    if ( GlobalConfig.MHD ) list ~= ["B.x", "B.y", "B.z"];
    list ~= ["p", "a", "mu"];
    auto gm = GlobalConfig.gmodel;
    foreach(i; 0 .. gm.n_modes) list ~= "k[" ~ to!string(i) ~ "]";
    list ~= ["mu_t", "k_t", "S"];
    if ( GlobalConfig.radiation ) list ~= ["Q_rad_org", "f_rad_org", "Q_rE_rad"];
    list ~= ["tke", "omega"];
    foreach(i; 0 .. gm.n_species) {
	auto name = cast(char[]) gm.species_name(i);
	// Clean up name, replacing internal spaces with dashes.
	name = strip(name);
	auto indx = indexOf(name, ' ');
	while ( indx > 0 ) { name[indx] = '-'; indx = indexOf(name, ' '); } 
	list ~= ["massf[" ~ to!string(i) ~ "]-" ~ to!string(name)];
    }
    if ( gm.n_species > 1 ) list ~= ["dt_chem"];
    foreach(i; 0 .. gm.n_modes) list ~= ["e[" ~ to!string(i) ~ "]", "T[" ~ to!string(i) ~ "]"];
    if ( gm.n_modes > 1 ) list ~= ["dt_therm"];
    return list;
} // end variable_list_for_cell()
